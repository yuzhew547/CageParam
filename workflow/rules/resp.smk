# resp.smk - optional RESP charge derivation via Multiwfn.
#
# Only fires when the user invoked QM separately and parked the resulting
# molden at scr_1/inputfile.molden but did NOT pre-place inputfile.chg.
# Heavy step (hours on a 16-core node). Multiwfn must be installed separately
# (closed-source), with its path in $CAGEPIPE_MULTIWFN or via --config multiwfn=.

rule respfit:
    input:
        molden = _join(CAGE, "scr_1/inputfile.molden"),
    output:
        chg = _join(CAGE, "inputfile.chg"),
    params:
        multiwfn = MULTIWFN,
        cage = CAGE or ".",
    threads: 16
    shell:
        r"""
        cd "{params.cage}"
        if [ -n "{params.multiwfn}" ]; then
            cagepipe-respfit scr_1/inputfile.molden -o inputfile.chg --keep-aux \
                --multiwfn "{params.multiwfn}" --nthreads {threads}
        else
            cagepipe-respfit scr_1/inputfile.molden -o inputfile.chg --keep-aux \
                --nthreads {threads}
        fi
        """
