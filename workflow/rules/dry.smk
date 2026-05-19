# dry.smk - core parametrization chain:
#   inputfile.xyz + inputfile.chg
#     -> pdb4munro     -> bone.pdb + per-residue mol2s
#     -> munro         -> munro.frcmod  (also runs antechamber+parmchk2 per template)
#     -> tleapgen      -> tleap.in
#     -> tleap         -> {prefix}_dry.{pdb,prmtop,inpcrd}

rule pdb4munro:
    """xyz + chg -> bone.pdb + LA*.mol2, P*.mol2, ... + Ltemp*_template.pdb"""
    input:
        xyz = _join(CAGE, "inputfile.xyz"),
        chg = _join(CAGE, "inputfile.chg"),
    output:
        bone = _join(CAGE, "bone.pdb"),
    params:
        cage = CAGE or ".",
    shell:
        r"""
        cd "{params.cage}"
        pdb4munro inputfile.xyz --chg inputfile.chg -o bone.pdb
        """


rule munro:
    """bone.pdb -> munro.frcmod (auto-runs antechamber+parmchk2 on templates).
    The GAFF data file defaults to the bundled gaff2.dat when --gaff is blank."""
    input:
        bone = _join(CAGE, "bone.pdb"),
    output:
        frcmod = _join(CAGE, "munro.frcmod"),
    params:
        gaff = GAFF,
        cage = CAGE or ".",
    shell:
        r"""
        cd "{params.cage}"
        if [ -n "{params.gaff}" ]; then
            munro -p bone.pdb -g "{params.gaff}" --auto-from-pdb -o munro.frcmod
        else
            munro -p bone.pdb --auto-from-pdb -o munro.frcmod
        fi
        """


rule tleapgen:
    """bone.pdb + munro.frcmod -> tleap.in"""
    input:
        bone   = _join(CAGE, "bone.pdb"),
        frcmod = _join(CAGE, "munro.frcmod"),
    output:
        leapin = _join(CAGE, "tleap.in"),
    params:
        prefix = PREFIX,
        cage = CAGE or ".",
    shell:
        r"""
        cd "{params.cage}"
        tleapgen -p bone.pdb -o tleap.in --prefix {params.prefix}
        """


rule tleap:
    """tleap.in -> {prefix}_dry.{pdb,prmtop,inpcrd}"""
    input:
        leapin = _join(CAGE, "tleap.in"),
        frcmod = _join(CAGE, "munro.frcmod"),
        bone   = _join(CAGE, "bone.pdb"),
    output:
        pdb    = _join(CAGE, f"{PREFIX}_dry.pdb"),
        prmtop = _join(CAGE, f"{PREFIX}_dry.prmtop"),
        inpcrd = _join(CAGE, f"{PREFIX}_dry.inpcrd"),
    params:
        cage = CAGE or ".",
    shell:
        r"""
        cd "{params.cage}"
        tleap -s -f tleap.in
        """
