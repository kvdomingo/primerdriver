from enum import Enum
from json import load
from typing import Any

from numpy import array
from pandas import DataFrame
from tabulate import tabulate

from primerdriver.checks import SequenceChecks
from primerdriver.config import BASE_DIR, get_settings
from primerdriver.log import logger

settings = get_settings()


class OperationMode(Enum):
    CHARACTERIZATION = "CHAR"
    DNA = "DNA"
    PROTEIN = "PRO"


class MutationType(Enum):
    SUBSTITUTION = "S"
    INSERTION = "I"
    DELETION = "D"


class PrimerMode(Enum):
    COMPLEMENTARY = "complementary"
    OVERLAPPING = "overlapping"


class PrimerDesign:
    df: DataFrame

    def __init__(
        self,
        mode: str,
        sequence,
        mutation_type,
        mismatched_bases=None,
        target=None,
        replacement=None,
        position=None,
        savename=None,
        print_buffer=20,
        config_file=BASE_DIR / "primerdriver" / "settings.json",  # deprecated
        **kwargs,
    ):
        self.settings: dict[str, Any] = settings
        self.mode = OperationMode(mode.upper())
        self.sequence: str = sequence.upper()
        self.mutation_type = MutationType(mutation_type[0].upper())
        self.mismatched_bases: int | None = (
            int(mismatched_bases) if mismatched_bases is not None else None
        )
        self.replacement: str | None = (
            replacement.upper() if replacement is not None else None
        )
        self.position: int | None = int(position) - 1 if position is not None else None
        self.target: str | None = target.upper() if target is not None else None
        self.Tm_range: tuple[float, float] = (
            self.settings["Tm_range_min"],
            self.settings["Tm_range_max"],
        )
        self.length_range: tuple[int, int] = (
            self.settings["length_min"],
            self.settings["length_max"],
        )
        self.gc_range: tuple[float, float] = (
            self.settings["gc_range_min"],
            self.settings["gc_range_max"],
        )
        self.flank5_range: tuple[int, int] = (
            self.settings["flank5_range_min"],
            self.settings["flank5_range_max"],
        )
        self.flank3_range: tuple[int, int] = (
            self.settings["flank3_range_min"],
            self.settings["flank3_range_max"],
        )
        self.forward_overlap5: int = self.settings["forward_overlap5"]
        self.forward_overlap3: int = self.settings["forward_overlap3"]
        self.terminate_gc = bool(self.settings["terminate_gc"])
        self.center_mutation = bool(self.settings["center_mutation"])
        self.primer_mode = PrimerMode(self.settings["primer_mode"])
        self.print_buffer: int = print_buffer
        self.expression_name: str = self.settings["expression_system"]
        self.savename = savename
        with open(BASE_DIR / "primerdriver" / "lut.json", "r", encoding="utf-8") as f:
            self.lut: dict[str, Any] = load(f)
        with open(
            BASE_DIR
            / "primerdriver"
            / "expression_system"
            / f"{self.expression_name}.json",
            "r",
        ) as f:
            self.expression_system: dict[str, str] = load(f)

    @staticmethod
    def calculate_gc_content(seq: str | list[str]) -> float:
        return (seq.count("G") + seq.count("C")) / len(seq)

    @staticmethod
    def calculate_mismatch(seq: str | list[str], mismatched_bases: int) -> float:
        return mismatched_bases / len(seq)

    def get_reverse_complement(self, seq: str | list[str]) -> list[str]:
        seq = list(seq)
        return [self.lut["complement"][b] for b in seq][::-1]

    @staticmethod
    def is_gc_end(sequence: str | list[str]) -> bool:
        sequence = "".join(sequence)
        return (sequence.startswith("G") or sequence.startswith("C")) and (
            sequence.endswith("G") or sequence.endswith("C")
        )

    @staticmethod
    def calculate_melting_temperature(
        seq: str | list[str],
        mutation_type: MutationType,
        replacement: str | list[str],
        gc_content: float,
        mismatch: float,
    ) -> float:
        gc_content = int(gc_content * 100)
        mismatch = int(mismatch * 100)
        if mutation_type == MutationType.SUBSTITUTION:
            sequence_length = len(seq)
            return 81.5 + 0.41 * gc_content - 675 / sequence_length - mismatch
        else:
            if mutation_type == MutationType.INSERTION:
                sequence_length = len(seq)
            elif mutation_type == MutationType.DELETION:
                sequence_length = len(seq)
            else:
                sequence_length = len(seq) - len(replacement)
            return 81.5 + 0.41 * gc_content - 675 / sequence_length

    def characterize_primer(
        self,
        sequence: str | list[str],
        mutation_type: MutationType,
        replacement: str | list[str] | None,
        mismatched_bases: int,
        index: int = None,
        reverse: list[str] = None,
    ) -> DataFrame:
        mol_weight = self.lut["mol_weight"]
        seq = list(sequence)
        if (
            self.primer_mode == PrimerMode.COMPLEMENTARY
            or self.mode == OperationMode.CHARACTERIZATION
        ):
            rev = "".join(self.get_reverse_complement(sequence))
        else:
            rev = reverse

        forwards = {
            "sequence": sequence,
            "length": len(seq),
            "gc_content": self.calculate_gc_content(seq),
            "mismatch": self.calculate_mismatch(seq, mismatched_bases),
            "gc_end": self.is_gc_end(seq),
            "mol_weight": sum(float(mol_weight[b]) for b in seq),
        }
        forwards["Tm"] = self.calculate_melting_temperature(
            seq,
            mutation_type,
            replacement,
            forwards["gc_content"],
            forwards["mismatch"],
        )

        reverses = {
            "sequence": rev,
            "length": len(rev),
            "gc_content": self.calculate_gc_content(rev),
            "mismatch": self.calculate_mismatch(rev, mismatched_bases),
            "gc_end": self.is_gc_end(rev),
            "mol_weight": sum(float(mol_weight[b]) for b in rev),
        }
        reverses["Tm"] = self.calculate_melting_temperature(
            seq,
            mutation_type,
            replacement,
            reverses["gc_content"],
            reverses["mismatch"],
        )

        col = [
            "Forward",
            "Reverse",
            "Fwd length",
            "Rev length",
            "Fwd GC content",
            "Rev GC content",
            "Fwd melting temp",
            "Rev melting temp",
            "Fwd mol. weight",
            "Rev mol. weight",
            "Fwd mismatch",
            "Rev mismatch",
            "Fwd ends in G/C",
            "Rev ends in G/C",
        ]
        dat = [
            forwards["sequence"],
            reverses["sequence"],
            f'{forwards["length"]} bp',
            f'{reverses["length"]} bp',
            f'{forwards["gc_content"]*100:.2f}%',
            f'{reverses["gc_content"]*100:.2f}%',
            f'{forwards["Tm"]:.2f} C',
            f'{reverses["Tm"]:.2f} C',
            f'{forwards["mol_weight"]:.2f} g/mol',
            f'{reverses["mol_weight"]:.2f} g/mol',
            f'{forwards["mismatch"]*100:.2f}%',
            f'{reverses["mismatch"]*100:.2f}%',
            forwards["gc_end"],
            reverses["gc_end"],
        ]
        if index is None:
            index = 1
        if index - 1 < self.print_buffer:
            print(
                "\n",
                tabulate(
                    array([col, dat]).T, headers=[f"Primer {index}"], tablefmt="orgtbl"
                ),
                sep="",
            )
        elif index - 1 == self.print_buffer:
            logger.info("Too many results; truncating output...")
        dat = array([dat])
        df = DataFrame(data=dat, columns=col, index=[index])
        self.df = df
        return df

    def substitution(  # noqa: C901
        self,
        sequence: str | list[str],
        mutation_type: MutationType,
        target: str | list[str] | None,
        replacement: str | list[str],
        start_position: int,
        mismatched_bases: int,
    ) -> DataFrame | None:
        if self.primer_mode == PrimerMode.COMPLEMENTARY:
            valid_primers = []
            seq = list(sequence)
            sequence_length = len(replacement)
            seq[start_position - 1] = replacement
            for f5 in range(*self.flank5_range):
                for f3 in range(*self.flank3_range):
                    if abs(f5 - f3) > 1 and self.center_mutation:
                        continue
                    candidate1 = seq[start_position - 1 - f5 : start_position - 1]
                    candidate3 = list(replacement)
                    candidate2 = seq[
                        start_position
                        + sequence_length
                        - 1 : start_position
                        + sequence_length
                        + f3
                    ]
                    candidate = candidate1 + candidate3 + candidate2
                    candidate = "".join(candidate)
                    if len(candidate) == 0:
                        continue
                    gc_content = self.calculate_gc_content(candidate)
                    mismatch = self.calculate_mismatch(candidate, mismatched_bases)
                    melting_temperature = self.calculate_melting_temperature(
                        candidate, mutation_type, replacement, gc_content, mismatch
                    )
                    sc = SequenceChecks(candidate)
                    valid_gc = sc.check_gc_content(self.gc_range)
                    valid_temp = sc.check_Tm(melting_temperature, self.Tm_range)
                    valid_ends = sc.check_ends_gc(self.terminate_gc)
                    valid_length = sc.check_sequence_length(self.length_range)
                    if valid_gc and valid_temp and valid_ends and valid_length:
                        valid_primers.append(candidate)
        else:
            valid_primers = {}
            forward_sequence = list(sequence)
            sequence_length = len(replacement)
            forward_sequence[start_position - 1] = replacement
            for f5 in range(*self.flank5_range):
                for f3 in range(*self.flank3_range):
                    if abs(f5 - f3) > 1 and self.center_mutation:
                        continue
                    candidate1 = forward_sequence[
                        start_position - 1 - f5 : start_position - 1
                    ]
                    candidate3 = list(replacement)
                    candidate2 = forward_sequence[
                        start_position
                        + sequence_length
                        - 1 : start_position
                        + sequence_length
                        + f3
                    ]
                    candidate = candidate1 + candidate3 + candidate2
                    candidate = "".join(candidate)
                    if len(candidate) == 0:
                        continue
                    gc_content = self.calculate_gc_content(candidate)
                    mismatch = self.calculate_mismatch(candidate, mismatched_bases)
                    melting_temperature = self.calculate_melting_temperature(
                        candidate, mutation_type, replacement, gc_content, mismatch
                    )
                    sc = SequenceChecks(candidate)
                    valid_gc = sc.check_gc_content(self.gc_range)
                    valid_temp = sc.check_Tm(melting_temperature, self.Tm_range)
                    valid_ends = sc.check_ends_gc(self.terminate_gc)
                    valid_length = sc.check_sequence_length(self.length_range)
                    if valid_gc and valid_temp and valid_ends and valid_length:
                        valid_primers[candidate] = []
            for primers in valid_primers.keys():
                sequence = (
                    sequence[: start_position - 1]
                    + replacement
                    + sequence[start_position - 1 + sequence_length :]
                )
                start = sequence.find(primers)
                end = start + len(primers) - 1
                while (
                    start < self.position - self.forward_overlap5
                    and end > self.position + sequence_length + self.forward_overlap3
                ):
                    start -= 1
                    end -= 1
                for i in range(self.position - self.forward_overlap5, end):
                    for j in range(self.flank3_range[1]):
                        candidate = sequence[start - j : end]
                        if len(candidate) == 0:
                            continue
                        gc_content = self.calculate_gc_content(candidate)
                        mismatch = self.calculate_mismatch(candidate, mismatched_bases)
                        forward_melting_temp = melting_temperature
                        melting_temperature = self.calculate_melting_temperature(
                            candidate, mutation_type, replacement, gc_content, mismatch
                        )
                        sc = SequenceChecks(candidate)
                        valid_gc = sc.check_gc_content(self.gc_range)
                        valid_temp = sc.check_Tm(melting_temperature, self.Tm_range)
                        valid_trange = sc.check_close_Tm(
                            forward_melting_temp, melting_temperature
                        )
                        valid_ends = sc.check_ends_gc(self.terminate_gc)
                        valid_length = sc.check_sequence_length(self.length_range)
                        if (
                            valid_gc
                            and valid_temp
                            and valid_ends
                            and valid_length
                            and valid_trange
                        ):
                            valid_primers[primers].append(candidate)
                    end += 1

        if len(valid_primers) == 0 or (
            isinstance(valid_primers, dict)
            and len([val for values in list(valid_primers.values()) for val in values])
            == 0
        ):
            print("No valid primers found")
            return
        else:
            df = []
            print(f"\nGenerated forward primers: {len(valid_primers)}")
            if self.mode == OperationMode.PROTEIN:
                print(f"Using expression system: {self.expression_name}")
            if self.primer_mode == PrimerMode.COMPLEMENTARY:
                for i, p in enumerate(valid_primers):
                    df.append(
                        self.characterize_primer(
                            p, mutation_type, replacement, mismatched_bases, i + 1
                        )
                    )
            else:
                count = 1
                for _, (k, v) in enumerate(valid_primers.items()):
                    for r in v:
                        df.append(
                            self.characterize_primer(
                                k,
                                mutation_type,
                                replacement,
                                mismatched_bases,
                                count,
                                r,
                            )
                        )
                        count += 1
        return df

    def deletion(  # noqa: C901
        self,
        sequence: str | list[str],
        mutation_type: MutationType,
        target: str | list[str] | None,
        replacement: str | list[str],
        start_position: int,
        mismatched_bases: int,
    ) -> DataFrame | None:
        if self.primer_mode == PrimerMode.COMPLEMENTARY:
            valid_primers = []
            sequence_length = len(self.target)
            seq = list(sequence)
            for f5 in range(*self.flank5_range):
                for f3 in range(*self.flank3_range):
                    if abs(f5 - f3) > 1 and self.center_mutation:
                        continue
                    if self.mode == OperationMode.DNA:
                        candidate1 = seq[start_position - 1 - f5 : start_position - 1]
                        candidate2 = seq[
                            start_position
                            + sequence_length
                            - 1 : start_position
                            + sequence_length
                            + f3
                        ]
                    elif self.mode == OperationMode.PROTEIN:
                        candidate1 = seq[start_position - 1 - f5 : start_position - 1]
                        candidate2 = seq[
                            start_position
                            + sequence_length * 3
                            - 1 : start_position
                            + sequence_length
                            + f3
                        ]
                    candidate = candidate1 + candidate2
                    candidate = "".join(candidate)
                    if len(candidate) == 0:
                        continue
                    gc_content = self.calculate_gc_content(candidate)
                    mismatch = self.calculate_mismatch(candidate, mismatched_bases)
                    melting_temperature = self.calculate_melting_temperature(
                        candidate, mutation_type, replacement, gc_content, mismatch
                    )
                    sc = SequenceChecks(candidate)
                    valid_gc = sc.check_gc_content(self.gc_range)
                    valid_temp = sc.check_Tm(melting_temperature, self.Tm_range)
                    valid_ends = sc.check_ends_gc(self.terminate_gc)
                    valid_length = sc.check_sequence_length(self.length_range)
                    if valid_gc and valid_temp and valid_ends and valid_length:
                        valid_primers.append(candidate)
        else:
            valid_primers = {}
            sequence_length = len(self.target)
            seq = list(sequence)
            for f5 in range(*self.flank5_range):
                for f3 in range(*self.flank3_range):
                    if abs(f5 - f3) > 1 and self.center_mutation:
                        continue
                    if self.mode == OperationMode.DNA:
                        candidate1 = seq[start_position - 1 - f5 : start_position - 1]
                        candidate2 = seq[
                            start_position
                            + sequence_length
                            - 1 : start_position
                            + sequence_length
                            + f3
                        ]
                    else:
                        candidate1 = seq[start_position - 1 - f5 : start_position - 1]
                        candidate2 = seq[
                            start_position
                            + sequence_length * 3
                            - 1 : start_position
                            + sequence_length
                            + f3
                        ]
                    candidate = candidate1 + candidate2
                    candidate = "".join(candidate)
                    if len(candidate) == 0:
                        continue
                    gc_content = self.calculate_gc_content(candidate)
                    mismatch = self.calculate_mismatch(candidate, mismatched_bases)
                    melting_temperature = self.calculate_melting_temperature(
                        candidate, mutation_type, replacement, gc_content, mismatch
                    )
                    sc = SequenceChecks(candidate)
                    valid_gc = sc.check_gc_content(self.gc_range)
                    valid_temp = sc.check_Tm(melting_temperature, self.Tm_range)
                    valid_ends = sc.check_ends_gc(self.terminate_gc)
                    valid_length = sc.check_sequence_length(self.length_range)
                    if valid_gc and valid_temp and valid_ends and valid_length:
                        valid_primers[candidate] = []
            if self.mode == OperationMode.DNA:
                sequence = (
                    sequence[: start_position - 1]
                    + sequence[start_position + sequence_length - 1 :]
                )
            elif self.mode == OperationMode.PROTEIN:
                sequence = (
                    sequence[: start_position - 1]
                    + sequence[start_position + sequence_length * 3 - 1 :]
                )
            for primers in valid_primers:
                start = sequence.find(primers)
                end = start + len(primers) - 1
                while (
                    start < self.position - self.forward_overlap5
                    and end > self.position + sequence_length + self.forward_overlap3
                ):
                    start -= 1
                    end -= 1
                for i in range(self.position - self.forward_overlap5, end):
                    for j in range(self.flank3_range[1]):
                        candidate = sequence[start - j : end]
                        if len(candidate) == 0:
                            continue
                        gc_content = self.calculate_gc_content(candidate)
                        mismatch = self.calculate_mismatch(candidate, mismatched_bases)
                        forward_melting_temp = melting_temperature
                        melting_temperature = self.calculate_melting_temperature(
                            candidate, mutation_type, replacement, gc_content, mismatch
                        )
                        sc = SequenceChecks(candidate)
                        valid_gc = sc.check_gc_content(self.gc_range)
                        valid_temp = sc.check_Tm(melting_temperature, self.Tm_range)
                        valid_trange = sc.check_close_Tm(
                            forward_melting_temp, melting_temperature
                        )
                        valid_ends = sc.check_ends_gc(self.terminate_gc)
                        valid_length = sc.check_sequence_length(self.length_range)
                        if (
                            valid_gc
                            and valid_temp
                            and valid_ends
                            and valid_length
                            and valid_trange
                        ):
                            valid_primers[primers].append(candidate)
                    end += 1

        if not len(valid_primers) > 0:
            print("No valid primers found")
            return
        else:
            df = []
            print(valid_primers)
            print(f"\nGenerated forward primers: {len(valid_primers)}")
            if self.primer_mode == "complementary":
                for i, p in enumerate(valid_primers):
                    df.append(
                        self.characterize_primer(
                            p, mutation_type, replacement, mismatched_bases, i + 1
                        )
                    )
            else:
                count = 1
                for _, (k, v) in enumerate(valid_primers.items()):
                    for r in v:
                        df.append(
                            self.characterize_primer(
                                k,
                                mutation_type,
                                replacement,
                                mismatched_bases,
                                count,
                                r,
                            )
                        )
                        count += 1
        return df

    def insertion(  # noqa: C901
        self,
        sequence: str | list[str],
        mutation_type: MutationType,
        target: str | list[str] | None,
        replacement: str | list[str],
        start_position: int,
        mismatched_bases: int,
    ) -> DataFrame | None:
        if self.primer_mode == PrimerMode.COMPLEMENTARY:
            valid_primers = []
            seq = list(sequence)
            sequence_length = len(replacement)
            seq[start_position - 1] = replacement + seq[start_position - 1]
            for f5 in range(*self.flank5_range):
                for f3 in range(*self.flank3_range):
                    if abs(f5 - f3) > 1 and self.center_mutation:
                        continue
                    candidate1 = seq[start_position - 1 - f5 : start_position - 1]
                    candidate2 = seq[
                        start_position - 1 : start_position + sequence_length + f3
                    ]
                    candidate = candidate1 + candidate2
                    candidate = "".join(candidate)
                    if len(candidate) == 0:
                        continue
                    gc_content = self.calculate_gc_content(candidate)
                    mismatch = self.calculate_mismatch(candidate, mismatched_bases)
                    melting_temp = self.calculate_melting_temperature(
                        candidate, mutation_type, replacement, gc_content, mismatch
                    )
                    sc = SequenceChecks(candidate)
                    valid_gc = sc.check_gc_content(self.gc_range)
                    valid_temp = sc.check_Tm(melting_temp, self.Tm_range)
                    valid_ends = sc.check_ends_gc(self.terminate_gc)
                    valid_length = sc.check_sequence_length(self.length_range)
                    if valid_gc and valid_temp and valid_ends and valid_length:
                        valid_primers.append(candidate)
        else:
            valid_primers = dict()
            seq = list(sequence)
            sequence_length = len(replacement)
            seq[start_position - 1] = replacement + seq[start_position - 1]
            for f5 in range(*self.flank5_range):
                for f3 in range(*self.flank3_range):
                    if abs(f5 - f3) > 1 and self.center_mutation:
                        continue
                    candidate1 = seq[start_position - 1 - f5 : start_position - 1]
                    candidate2 = seq[
                        start_position - 1 : start_position + sequence_length + f3
                    ]
                    candidate = candidate1 + candidate2
                    candidate = "".join(candidate)
                    if len(candidate) == 0:
                        continue
                    gc_content = self.calculate_gc_content(candidate)
                    mismatch = self.calculate_mismatch(candidate, mismatched_bases)
                    melting_temp = self.calculate_melting_temperature(
                        candidate, mutation_type, replacement, gc_content, mismatch
                    )
                    sc = SequenceChecks(candidate)
                    valid_gc = sc.check_gc_content(self.gc_range)
                    valid_temp = sc.check_Tm(melting_temp, self.Tm_range)
                    valid_ends = sc.check_ends_gc(self.terminate_gc)
                    valid_length = sc.check_sequence_length(self.length_range)
                    if valid_gc and valid_temp and valid_ends and valid_length:
                        valid_primers[candidate] = []
            sequence = (
                sequence[: start_position - 1]
                + replacement
                + sequence[start_position - 1 :]
            )
            for primers in valid_primers:
                start = sequence.find(primers)
                end = start + len(primers) - 1
                while (
                    start < self.position - self.forward_overlap5
                    and end > self.position + sequence_length + self.forward_overlap3
                ):
                    start -= 1
                    end -= 1
                for i in range(self.position - self.forward_overlap5, end):
                    for j in range(self.flank3_range[1]):
                        candidate = sequence[start - j : end]
                        if len(candidate) == 0:
                            continue
                        gc_content = self.calculate_gc_content(candidate)
                        mismatch = self.calculate_mismatch(candidate, mismatched_bases)
                        fwd_Tm = melting_temp
                        melting_temp = self.calculate_melting_temperature(
                            candidate, mutation_type, replacement, gc_content, mismatch
                        )
                        sc = SequenceChecks(candidate)
                        valid_gc = sc.check_gc_content(self.gc_range)
                        valid_temp = sc.check_Tm(melting_temp, self.Tm_range)
                        valid_trange = sc.check_close_Tm(fwd_Tm, melting_temp)
                        valid_ends = sc.check_ends_gc(self.terminate_gc)
                        valid_length = sc.check_sequence_length(self.length_range)
                        if (
                            valid_gc
                            and valid_temp
                            and valid_ends
                            and valid_length
                            and valid_trange
                        ):
                            valid_primers[primers].append(candidate)
                    end += 1

        if not len(valid_primers) > 0:
            print("No valid primers found")
            return
        else:
            df = []
            print(f"\nGenerated forward primers: {len(valid_primers)}")
            if self.primer_mode == PrimerMode.COMPLEMENTARY:
                for i, p in enumerate(valid_primers):
                    df.append(
                        self.characterize_primer(
                            p, mutation_type, replacement, mismatched_bases, i + 1
                        )
                    )
            else:
                count = 1
                for _, (k, v) in enumerate(valid_primers.items()):
                    for r in v:
                        df.append(
                            self.characterize_primer(
                                k,
                                mutation_type,
                                replacement,
                                mismatched_bases,
                                count,
                                r,
                            )
                        )
                        count += 1
        return df

    def dna_based(self):
        if self.mutation_type == MutationType.DELETION:
            if self.mismatched_bases is None:
                self.mismatched_bases = len(self.target)
            result = self.deletion(
                self.sequence,
                self.mutation_type,
                self.target,
                self.replacement,
                self.position,
                self.mismatched_bases,
            )
        elif self.mutation_type == MutationType.INSERTION:
            if self.mismatched_bases is None:
                self.mismatched_bases = len(self.replacement)
            result = self.insertion(
                self.sequence,
                self.mutation_type,
                self.target,
                self.replacement,
                self.position,
                self.mismatched_bases,
            )
        else:
            if self.mismatched_bases is None:
                self.mismatched_bases = len(self.replacement)
            result = self.substitution(
                self.sequence,
                self.mutation_type,
                self.target,
                self.replacement,
                self.position,
                self.mismatched_bases,
            )
        return result

    def protein_based(self):
        if self.mutation_type == MutationType.DELETION:
            if self.mismatched_bases is None:
                self.mismatched_bases = len(self.target)
            rna = "".join([self.expression_system[b] for b in self.sequence])
            dna = rna.replace("U", "T")
            result = self.deletion(
                dna,
                self.mutation_type,
                self.target,
                self.replacement,
                self.position * 3,
                self.mismatched_bases,
            )
        elif self.mutation_type == MutationType.INSERTION:
            if self.mismatched_bases is None:
                self.mismatched_bases = len(self.replacement)
            rna = "".join([self.expression_system[b] for b in self.sequence])
            dna = rna.replace("U", "T")
            replacement = "".join(self.expression_system[self.replacement]).replace(
                "U", "T"
            )
            result = self.insertion(
                dna,
                self.mutation_type,
                self.target,
                replacement,
                self.position * 3,
                self.mismatched_bases,
            )
        else:
            if self.mismatched_bases is None:
                self.mismatched_bases = len(self.replacement)
            rna = "".join([self.expression_system[b] for b in self.sequence])
            dna = rna.replace("U", "T")
            target = "".join(self.expression_system[self.target]).replace("U", "T")
            replacement = "".join(self.expression_system[self.replacement]).replace(
                "U", "T"
            )
            result = self.substitution(
                dna,
                self.mutation_type,
                target,
                replacement,
                self.position * 3,
                self.mismatched_bases,
            )
        return result

    def main(self):
        if self.mode == OperationMode.DNA:
            df = self.dna_based()
        elif self.mode == OperationMode.PROTEIN:
            df = self.protein_based()
        else:
            df = self.characterize_primer(
                self.sequence,
                self.mutation_type,
                self.replacement,
                self.mismatched_bases,
            )
        self.df = df
