"""Beam search implementation for template-free models with KV caching.

Uses the TransformerDecoder's built-in KV cache infrastructure to avoid
re-decoding the full token sequence at every step. This reduces decoder work
from O(beam_size * N^2) to O(beam_size * N).
"""

import numpy as np
import torch
import torch.nn.functional as F
from dataclasses import dataclass
from typing import List


@dataclass
class BeamSearchResult:
    """A single beam search result."""
    token_ids: List[int]
    score: float  # log probability


def beam_search_decode(model, batch, vocab, beam_size=10, n_best=10,
                       max_length=512, temperature=1.0, device=None):
    """Run beam search decoding with KV caching.

    Args:
        model: Model with encode_and_reshape, decoder, output_layer.
        batch: Input batch.
        vocab: Dictionary mapping token strings to indices.
        beam_size: Number of beams.
        n_best: Number of best results to return.
        max_length: Maximum decoding length.
        temperature: Temperature for softmax.
        device: torch device.

    Returns:
        List of BeamSearchResult, sorted by score (best first).
    """
    sos_id = vocab["_SOS"]
    eos_id = vocab["_EOS"]

    with torch.no_grad():
        padded_memory_bank, memory_lengths = model.encode_and_reshape(batch)

        if beam_size == 1:
            return _greedy_decode(model, padded_memory_bank, memory_lengths,
                                  vocab, max_length, temperature, device)

        mem_len = padded_memory_bank.size(0)

        # Step 0: Process SOS token with batch=1 to initialize cache
        model.decoder.state["src"] = torch.zeros(mem_len, 1, device=device)
        model.decoder.state["cache"] = None

        sos_input = torch.tensor([[[sos_id]]], device=device)  # (1, 1, 1)
        dec_out, _ = model.decoder(
            tgt=sos_input,
            memory_bank=padded_memory_bank,
            memory_lengths=memory_lengths,
            step=0,
            mode="train")

        # Get logits from step 0
        last_out = model.output_layer(dec_out[-1:]) / temperature  # (1, 1, vocab)
        log_probs = F.log_softmax(last_out.squeeze(0).squeeze(0), dim=-1)

        # Select top beam_size tokens
        topk_scores, topk_ids = log_probs.topk(beam_size)  # (beam_size,)

        # Track token history: (steps_so_far, beam_size)
        # Start with SOS + first predicted token for each beam
        token_history = torch.full((2, beam_size), sos_id, dtype=torch.long, device=device)
        token_history[1] = topk_ids

        beam_scores = topk_scores  # (beam_size,)

        # Expand memory_bank from (mem, 1, H) to (mem, beam_size, H)
        padded_memory_bank = padded_memory_bank.repeat(1, beam_size, 1)
        memory_lengths = memory_lengths.repeat(beam_size)

        # Expand cache from batch=1 to batch=beam_size
        model.decoder.state["src"] = torch.zeros(mem_len, beam_size, device=device)
        model.decoder.map_state(
            lambda t, dim: t.repeat_interleave(beam_size, dim=dim))

        finished_beams = []
        # Track which beams are finished (score set to -inf to prevent expansion)
        is_finished = torch.zeros(beam_size, dtype=torch.bool, device=device)

        # Check if any step-0 beams already produced EOS
        eos_mask = topk_ids.eq(eos_id)
        if eos_mask.any():
            for idx in eos_mask.nonzero(as_tuple=False).squeeze(-1):
                i = idx.item()
                toks = token_history[:, i].tolist()
                finished_beams.append(BeamSearchResult(
                    token_ids=toks,
                    score=beam_scores[i].item() / len(toks)))
                is_finished[i] = True
            beam_scores[is_finished] = -1e9

        for step in range(1, max_length):
            if len(finished_beams) >= n_best:
                break
            if is_finished.all():
                break

            # Input: last token only, shape (1, beam_size, 1)
            last_tokens = token_history[-1].unsqueeze(0).unsqueeze(-1)  # (1, beam_size, 1)

            dec_out, _ = model.decoder(
                tgt=last_tokens,
                memory_bank=padded_memory_bank,
                memory_lengths=memory_lengths,
                step=step,
                mode="train")

            # dec_out: (1, beam_size, hidden) -> logits
            last_out = model.output_layer(dec_out) / temperature  # (1, beam_size, vocab)
            log_probs = F.log_softmax(last_out.squeeze(0), dim=-1)  # (beam_size, vocab)
            vocab_size = log_probs.size(-1)

            # New scores: (beam_size, vocab)
            new_scores = beam_scores.unsqueeze(1) + log_probs

            # Flatten and select top beam_size
            flat_scores = new_scores.view(-1)  # (beam_size * vocab)
            topk_scores, topk_flat_ids = flat_scores.topk(beam_size)

            parent_indices = topk_flat_ids // vocab_size  # which beam each came from
            token_indices = topk_flat_ids % vocab_size  # which token was selected

            # Reorder cache to match new beam parents
            model.decoder.map_state(
                lambda t, dim: t.index_select(dim, parent_indices))

            # Reorder memory_bank
            padded_memory_bank = padded_memory_bank.index_select(1, parent_indices)

            # Reorder token history and append new tokens
            token_history = token_history[:, parent_indices]  # reorder columns
            token_history = torch.cat([
                token_history,
                token_indices.unsqueeze(0)
            ], dim=0)

            beam_scores = topk_scores
            is_finished = is_finished[parent_indices]

            # Check for EOS
            eos_mask = token_indices.eq(eos_id)
            if eos_mask.any():
                for idx in eos_mask.nonzero(as_tuple=False).squeeze(-1):
                    i = idx.item()
                    if not is_finished[i]:
                        toks = token_history[:, i].tolist()
                        finished_beams.append(BeamSearchResult(
                            token_ids=toks,
                            score=beam_scores[i].item() / len(toks)))
                        is_finished[i] = True
                beam_scores[is_finished] = -1e9

        # Add remaining alive beams
        for i in range(beam_size):
            if not is_finished[i]:
                toks = token_history[:, i].tolist()
                finished_beams.append(BeamSearchResult(
                    token_ids=toks,
                    score=beam_scores[i].item() / len(toks)))

        finished_beams.sort(key=lambda x: x.score, reverse=True)
        return finished_beams[:n_best]


def _greedy_decode(model, padded_memory_bank, memory_lengths,
                   vocab, max_length, temperature, device):
    """Greedy decoding with KV caching."""
    sos_id = vocab["_SOS"]
    eos_id = vocab["_EOS"]

    mem_len = padded_memory_bank.size(0)
    model.decoder.state["src"] = torch.zeros(mem_len, 1, device=device)
    model.decoder.state["cache"] = None

    tokens = [sos_id]

    # Step 0: SOS token, initialize cache
    sos_input = torch.tensor([[[sos_id]]], device=device)  # (1, 1, 1)
    dec_out, _ = model.decoder(
        tgt=sos_input,
        memory_bank=padded_memory_bank,
        memory_lengths=memory_lengths,
        step=0,
        mode="train")

    last_out = model.output_layer(dec_out[-1:]) / temperature
    log_probs = F.log_softmax(last_out.squeeze(0).squeeze(0), dim=-1)
    next_token = log_probs.argmax().item()
    tokens.append(next_token)

    if next_token == eos_id:
        return [BeamSearchResult(token_ids=tokens, score=0.0)]

    # Steps 1+: feed only last token with step=N
    for step in range(1, max_length):
        dec_in = torch.tensor([[[next_token]]], device=device)  # (1, 1, 1)

        dec_out, _ = model.decoder(
            tgt=dec_in,
            memory_bank=padded_memory_bank,
            memory_lengths=memory_lengths,
            step=step,
            mode="train")

        last_out = model.output_layer(dec_out[-1:]) / temperature
        log_probs = F.log_softmax(last_out.squeeze(0).squeeze(0), dim=-1)
        next_token = log_probs.argmax().item()
        tokens.append(next_token)

        if next_token == eos_id:
            break

    return [BeamSearchResult(token_ids=tokens, score=0.0)]


# ============================================================
# Naive implementations (kept for equivalence testing)
# ============================================================

def _beam_search_decode_naive(model, batch, vocab, beam_size=10, n_best=10,
                              max_length=512, temperature=1.0, device=None):
    """Original naive beam search that re-decodes full sequence each step."""
    pad_id = vocab["_PAD"]
    sos_id = vocab["_SOS"]
    eos_id = vocab["_EOS"]

    with torch.no_grad():
        padded_memory_bank, memory_lengths = model.encode_and_reshape(batch)

        if beam_size == 1:
            return _greedy_decode_naive(model, padded_memory_bank, memory_lengths,
                                       vocab, max_length, temperature, device)

        model.decoder.state["src"] = np.zeros(padded_memory_bank.size(0))

        alive_beams = [(torch.tensor([sos_id], device=device), 0.0)]
        finished_beams = []

        for step in range(max_length):
            if not alive_beams:
                break

            all_candidates = []
            for beam_idx, (tokens, score) in enumerate(alive_beams):
                model.decoder.state["src"] = np.zeros(padded_memory_bank.size(0))
                model.decoder.state["cache"] = None

                dec_in = tokens.unsqueeze(1).unsqueeze(-1)

                dec_out, _ = model.decoder(
                    tgt=dec_in,
                    memory_bank=padded_memory_bank,
                    memory_lengths=memory_lengths,
                    mode="train")

                last_out = model.output_layer(dec_out[-1:])
                last_out = last_out / temperature
                log_probs = F.log_softmax(last_out.squeeze(0).squeeze(0), dim=-1)

                topk_log_probs, topk_ids = log_probs.topk(beam_size)

                for k in range(beam_size):
                    new_token = topk_ids[k]
                    new_score = score + topk_log_probs[k].item()
                    new_tokens = torch.cat([tokens, new_token.unsqueeze(0)])
                    all_candidates.append((new_tokens, new_score))

            if not all_candidates:
                break

            all_candidates.sort(key=lambda x: x[1], reverse=True)
            alive_beams = all_candidates[:beam_size]

            new_alive = []
            for tokens, score in alive_beams:
                if tokens[-1].item() == eos_id:
                    finished_beams.append(BeamSearchResult(
                        token_ids=tokens.tolist(),
                        score=score / len(tokens)))
                else:
                    new_alive.append((tokens, score))
            alive_beams = new_alive

            if len(finished_beams) >= n_best:
                break

        for tokens, score in alive_beams:
            finished_beams.append(BeamSearchResult(
                token_ids=tokens.tolist(),
                score=score / len(tokens)))

        finished_beams.sort(key=lambda x: x.score, reverse=True)
        return finished_beams[:n_best]


def _greedy_decode_naive(model, padded_memory_bank, memory_lengths,
                         vocab, max_length, temperature, device):
    """Original naive greedy decoding that re-decodes full sequence each step."""
    sos_id = vocab["_SOS"]
    eos_id = vocab["_EOS"]

    model.decoder.state["src"] = np.zeros(padded_memory_bank.size(0))
    model.decoder.state["cache"] = None

    tokens = [sos_id]

    for step in range(max_length):
        dec_in = torch.tensor(tokens, device=device).unsqueeze(1).unsqueeze(-1)

        model.decoder.state["src"] = np.zeros(padded_memory_bank.size(0))
        model.decoder.state["cache"] = None

        dec_out, _ = model.decoder(
            tgt=dec_in,
            memory_bank=padded_memory_bank,
            memory_lengths=memory_lengths,
            mode="train")

        last_out = model.output_layer(dec_out[-1:])
        last_out = last_out / temperature
        log_probs = F.log_softmax(last_out.squeeze(0).squeeze(0), dim=-1)

        next_token = log_probs.argmax().item()
        tokens.append(next_token)

        if next_token == eos_id:
            break

    return [BeamSearchResult(token_ids=tokens, score=0.0)]
