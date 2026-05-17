"""cagepipe — metal-organic cage parametrization pipeline.

Pipeline (per-cage, runs inside the cage directory):
    inputfile.xyz + inputfile.chg
      -> pdb4munro      (xyz/chg -> bone.pdb + LA*.mol2, P*.mol2, ...)
      -> munro          (bone.pdb + GAFF -> frcmod)
      -> tleapgen       (bone.pdb -> tleap.in)
      -> tleap          (-> ori_dry.{pdb,prmtop,inpcrd})

Other entry points:
    cagepipe-respfit          Multiwfn RESP (whole-cage or differential)
    cagepipe-chgass           per-residue mol2 assignment from pre-charged templates
    cagepipe-seasoning        place N anions around the cage
    cagepipe-seasoning-inside place 1 PFOA inside the cavity (large-linker cages)
"""
from importlib.resources import files as _files

__version__ = "0.1.0"


def data_path(filename: str) -> str:
    """Absolute path to a file shipped under ``cagepipe/data/``.

    Example:
        >>> from cagepipe import data_path
        >>> data_path("gaff2.dat")
        '/.../site-packages/cagepipe/data/gaff2.dat'
    """
    return str(_files("cagepipe").joinpath("data", filename))
