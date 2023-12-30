# RetroRules pre-parsed files release rr02 for RetroPath 3 -- implicit Hs

## Context

More info at https://retrorules.org.

Pre-parsed rule files can be used as an alternative to the full
RetroRules SQLite database. Such files are ready to use with 
RetroPath Suite tools such as RetroPath2.0
(https://www.myexperiment.org/workflows/4987.html) and RetroPath3.0
(https://github.com/brsynth/RetroPath3).


### Licenses

This work is licensed under the Creative Commons Attribution 4.0
International License. To view a copy of this license, visit
http://creativecommons.org/licenses/by/4.0/.

The current release of RetroRules is based on MNXref v3.0 for the
chemical structures, the reaction definitions and their respective
crosslinks to other databases.

Please consider that MNXref information is gathered from different
metabolic databases (MetaCyc, KeGG, Bigg, Rhea, ...) and each have
different license policies (especially concerning the reuse for
commercial purpose). As a consequence, the license related to each
individual reaction rule depends of the licensing agreements of the
reference database(s) chosen in MNXref to describe each reaction as
well as the chemicals involved in it. This information is available
at https://www.metanetx.org/mnxdoc/mnxref.html.

Sequence annotations is obtained by cross-referencing MNXref v3.0
with Rhea release 81 and UniProt release 2017_04.


### Archive content

Present archive contains one file of rules:

- `retrorules_rr02_flat_all.tsv`: the complete dataset providing all the non-stereo rules from diameters 2 to 16


## Structure of file

Pre-parsed RetroRules file are TSV (tabulated separated values) having the following structure:

- `# Rule ID`: Reaction rule ID
- `Legacy ID`: Reaction rule ID according to the old (retrorules release 01) naming schema
- `Reaction_ID`: ID of the template reaction
- `Diameter`: Diameter of a rule
- `Reaction order`: The number of chemical structure on the left part of a reaction rule. By design in RetroRules all
reaction rules are mono-component rules, hence value is always `1`
- `Rule_SMARTS`: Reaction SMARTS encoding the chemical changes
- `Substrate_ID`: ID of the template substrate
- `Substrate_SMILES`: SMILES of the template substrate
- `Product_IDs`: IDs of the template products (separated by dot `.`)
- `Product_SMILES`: SMILES of the template products (separated by dot `.`)
- `Rule_SMILES`: Reaction rule expressed in the SMILES formalism
- `Rule_SMARTS_lite`: Reaction rule expressed in the SMARTS formalism but without the AAM (those are not executable)
- `Score`: Not normalized penalty score associated to
- `Score_normalized`: Normalized score. It ranges between 0 and 1, 1 is best. See RetroPath3 paper for details
on how this score is computed
- `Reaction_EC_number`: EC number annotations associated to the template reaction. Multiple annotations are separated
by `;`, `NOEC` is used when no EC number associated
- `Reaction_direction`: Directionality of the reference reaction (based on MetaCyc annotations). Possible values
are `1`, `0` and `-1` corresponding respectively to a left-to-right, reversible and right-to-left direction
- `Rule_relative_direction`: Direction of the rule comparatively to the reference reaction. Possible values are  `1`
and  `-1`. Rule direction is `1`   for a rule generated in the same direction than the reference reaction,  while `-1`
is used for a rules generated in the opposite direction. Reference reactions are always considered in their
left-to-right direction
- `Rule_usage`: Preferred usage for a given rule. Possible values are `retro`, `forward` and `both`. The preferred
usage takes into account (i) the directionality of the reaction of reference and (ii) the direction of the rule
relatively to its reference reaction. See below for details about the usage of rules.


## Reaction rule usages

`Forward` rules have been generated from reference reactions that are unidirectional and describe the transformation
in the natural direction, i.e. a substrate in a reference reaction is considered as substrate in a forward rule.
Forward rules can be helpful to identify promiscuous products of (bio)chemical reactions, or to predict sensing
enabling metabolic pathways.

By complementation, `retrosynthetic` rules have been generated in the reverse direction of unidirectional reactions,
i. e. a product in a reference reaction is considered as a substrate in the retrosynthetic rule. Such retrosynthetic
rules allow one to take advantage of retrosynthetic algorithm in order to design metabolic pathways for the production
of chemical compounds.

Finally, rules being predicted from bidirectionnal reactions can be used for `both` usages, since each compound of such
reactions can be used as a substrate or a product depending of which direction is considered.


### Explicit vs. implicit hydrogens

Hyrdrogen can be optionally expressed in the SMARTS reaction. If it
is the case, Hs are said 'explicit', while in the other case Hs are
said 'implicit' to contrast.