from typing import Literal

from pydantic import BaseModel

from primerdriver.primer_design import OperationMode


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


class DesignRequest(BaseModel):
    mode: OperationMode
    sequence: str
    mutation_type: str
    mismatched_bases: int | None = None
    target: str | None = None
    replacement: str | None = None
    position: int | None = None
    settings: AdvancedSettings | None = None


class VersionResponse(BaseModel):
    program_version: str
    web_version: str


class ExpressionSystemsResponse(BaseModel):
    data: list[str]


class MessageResponse(BaseModel):
    data: str


type DesignResult = dict[int | str, dict[int | str, str | int | float | bool]]
