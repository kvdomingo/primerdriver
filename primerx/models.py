from typing import Literal

from pydantic import BaseModel


class AdvancedSettings(BaseModel):
    Tm_range_min: float
    Tm_range_max: float
    gc_range_min: float
    gc_range_max: float
    length_min: int
    length_max: int
    flank5_range_min: int
    flank5_range_max: int
    flank3_range_min: int
    flank3_range_max: int
    forward_overlap5: int
    forward_overlap3: int
    terminate_gc: bool
    center_mutation: bool
    primer_mode: Literal["overlapping", "complementary"]
    expression_system: str
