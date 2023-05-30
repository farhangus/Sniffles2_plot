from dataclasses import dataclass

@dataclass
class VcfVariables():
    DEL:list
    DEL_GENOTYPE :list
    DEL_AF :list
    INS:list
    INS_GENOTYPE :list
    INS_AF :list
    DUP:list
    DUP_GENOTYPE :list
    DUP_AF :list
    INV:list
    INV_GENOTYPE :list
    INV_AF:list
    BND:list
    BND_GENOTYPE :list
    BND_AF :list
    @classmethod
    def new(cls):
        return cls([],[],[],[],[],[],[],[],[],[],[],[],[],[],[])