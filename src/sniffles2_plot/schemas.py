from dataclasses import dataclass


@dataclass
class GenomeChartData:
    points: list
    legend: str

    def count(self, labels):
        return [self.points.count(lbl) for lbl in labels]


@dataclass
class VcfVariables:
    DEL: list
    DEL_GENOTYPE: list
    DEL_AF: list
    INS: list
    INS_GENOTYPE: list
    INS_AF: list
    DUP: list
    DUP_GENOTYPE: list
    DUP_AF: list
    INV: list
    INV_GENOTYPE: list
    INV_AF: list
    BND: list
    BND_GENOTYPE: list
    BND_AF: list

    PHASED_DEL: list
    PHASED_DEL_GENOTYPE: list
    PHASED_DEL_AF: list
    PHASED_INS: list
    PHASED_INS_GENOTYPE: list
    PHASED_INS_AF: list
    PHASED_DUP: list
    PHASED_DUP_GENOTYPE: list
    PHASED_DUP_AF: list
    PHASED_INV: list
    PHASED_INV_GENOTYPE: list
    PHASED_INV_AF: list
    PHASED_BND: list
    PHASED_BND_GENOTYPE: list
    PHASED_BND_AF: list
    HAS_PHASED: bool=False

    @classmethod
    def new(cls):
        return cls([], [], [], [], [], [], [], [], [], [], [], [], [], [], [],[], [], [], [], [], [], [], [], [], [], [], [], [], [], [])
