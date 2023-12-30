import time
from functools import lru_cache

from rbc2.reaction_evaluation.scscore.standalone_model_numpy import sc_scorer

@lru_cache(maxsize=10000)
def get_complexity(smi, round_to=4):
    (smile, score) = sc_scorer.get_score_from_smi(smi)
    return round(score, round_to)

if __name__ == '__main__':
    smis = ['O=C(OC(C)(C)C)C(C=C1)=CC=C1CN(CC2)CCC2CN[C@H]3[C@@H](C3)C4=CC=CC=C4',
            'C#C[C@]1([C@H](C[C@@H](O1)N2C=NC3=C(N=C(N=C32)F)N)O)CO'
            'N[C@@H](CC(N1CCN2C(C1)=NN=C2C(F)(F)F)=O)CC3=C(F)C=C(F)C(F)=C3',
            'O=C(N[C@@H]1CCN(C)[C@@H](C2=NC(C=CC=C3)=C3N2)C1)CC4=CC=C(C#N)C=C4',
            'O=C(OC(C)(C)C)C(C=C1)=CC=C1CN(CC2)CCC2CN[C@H]3[C@@H](C3)C4=CC=CC=C4']

    times = []
    for i in range(500):
        for smi in smis:
            t0 = time.time()
            complexity = get_complexity(smi)
            t1 = time.time()
            time_taken = round(t1-t0,3)
            times.append(time_taken)
            print(f'Complexity in {time_taken} seconds')

    average_time = round(sum(times) / len(times), 3)
    print(f"Average time = {average_time}")
