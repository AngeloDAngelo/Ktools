"""Generate fictional tutorial sequences; no study data are used."""
from pathlib import Path
import random


def generate(destination):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    rng = random.Random(42)
    for group, count in [("reference", 24), ("control", 48)]:
        rows = []
        for i in range(count):
            sequence = "".join(rng.choices("ACGT", k=600))
            if group == "reference":
                # Replace a variable central segment with an ACG repeat.
                segment = "ACG" * (18 + i % 7)
                sequence = sequence[:180] + segment + sequence[180 + len(segment):]
            rows.append(f">{group}_{i:03d}\n{sequence}\n")
        (destination / f"{group}.fa").write_text("".join(rows))


if __name__ == "__main__":
    generate(Path(__file__).resolve().parent / "data")
