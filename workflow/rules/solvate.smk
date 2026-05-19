# solvate.smk - optional guest staging.
#
# These rules are NOT wired into `rule all` by default. Invoke explicitly:
#   snakemake --cores 1 seasoning_outside    # bone.pdb -> tastybone.pdb (N anions surround)
#   snakemake --cores 1 seasoning_inside     # bone.pdb -> tastybone.pdb (1 anion inside)
#
# The number of anions, anion template, etc. are config-driven.

ANION_TEMPLATE   = config.get("anion_template", "BFA")     # base name; needs <name>.{pdb,mol2,frcmod}
N_ANIONS_OUTSIDE = config.get("n_anions_outside", 24)


rule seasoning_outside:
    """Surround the cage with N anion copies (random shells outside the cavity)."""
    input:
        bone     = _join(CAGE, "bone.pdb"),
        tmpl_pdb = _join(CAGE, f"{ANION_TEMPLATE}.pdb"),
        tmpl_mol = _join(CAGE, f"{ANION_TEMPLATE}.mol2"),
        tmpl_frc = _join(CAGE, f"{ANION_TEMPLATE}.frcmod"),
    output:
        tasty = _join(CAGE, "tastybone.pdb"),
    params:
        n = N_ANIONS_OUTSIDE,
        cage = CAGE or ".",
        templ = ANION_TEMPLATE,
    shell:
        r"""
        cd "{params.cage}"
        # seasoning.py is interactive; feed it sensible defaults via heredoc.
        seasoning <<EOF
bone.pdb
tastybone.pdb
{params.n}
{params.templ}.pdb
{params.templ}.frcmod
{params.templ}.mol2
EOF
        """


rule seasoning_inside:
    """Place one anion (PFOA template) inside the cage cavity."""
    input:
        bone     = _join(CAGE, "bone.pdb"),
        tmpl_pdb = _join(CAGE, f"{ANION_TEMPLATE}.pdb"),
        tmpl_mol = _join(CAGE, f"{ANION_TEMPLATE}.mol2"),
        tmpl_frc = _join(CAGE, f"{ANION_TEMPLATE}.frcmod"),
    output:
        tasty = _join(CAGE, "tastybone_inside.pdb"),
    params:
        cage = CAGE or ".",
        templ = ANION_TEMPLATE,
    shell:
        r"""
        cd "{params.cage}"
        filling --cage bone.pdb --out tastybone_inside.pdb \
            --templ {params.templ}.pdb --resname {params.templ}
        """
