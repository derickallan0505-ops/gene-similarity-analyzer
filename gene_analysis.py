from itertools import combinations
from typing import Dict, List, Tuple

def read_fasta(file_path: str) -> Dict[str, str]:
    """Lê um arquivo FASTA e retorna um dicionário: {nome_do_organismo: sequência}."""
    sequences = {}
    header = None
    current = []

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sequences[header] = "".join(current).upper()
                header = line[1:].strip()
                current = []
            else:
                current.append(line)

    if header is not None:
        sequences[header] = "".join(current).upper()

    return sequences


def hamming_similarity(seq1: str, seq2: str) -> float:
    """Retorna identidade percentual por comparação direta (mesmo tamanho)."""
    if len(seq1) != len(seq2):
        raise ValueError("Sequências com tamanhos diferentes para comparação direta.")

    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return (matches / len(seq1)) * 100 if seq1 else 0.0


def global_identity(seq1: str, seq2: str) -> float:
    """
    Calcula similaridade percentual usando alinhamento global simples.
    Para genes de tamanhos parecidos, funciona bem.
    """
    if not seq1 or not seq2:
        return 0.0

    # Alinhamento simples: compara posições correspondentes até o menor comprimento
    min_len = min(len(seq1), len(seq2))
    matches = 0

    for i in range(min_len):
        if seq1[i] == seq2[i]:
            matches += 1

    identity = (matches / min_len) * 100
    return identity


def pairwise_similarity(sequences: Dict[str, str]) -> List[Tuple[str, str, float]]:
    """Calcula a similaridade entre todas as combinações de organismos."""
    result = []
    items = list(sequences.items())

    for (org1, seq1), (org2, seq2) in combinations(items, 2):
        identity = global_identity(seq1, seq2)
        result.append((org1, org2, round(identity, 2)))

    return result


def similarity_matrix(sequences: Dict[str, str]) -> Dict[str, Dict[str, float]]:
    """Cria uma matriz de similaridade entre todos os organismos."""
    matrix = {org: {} for org in sequences}
    items = list(sequences.items())

    for i, (org1, seq1) in enumerate(items):
        for org2, seq2 in items[i + 1:]:
            sim = global_identity(seq1, seq2)
            matrix[org1][org2] = round(sim, 2)
            matrix[org2][org1] = round(sim, 2)

    for org in sequences:
        matrix[org][org] = 100.0

    return matrix
