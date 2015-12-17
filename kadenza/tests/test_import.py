"""Basic test: can we import kadenza?"""

def test_import():
    from .. import kadenza
    try:
        kadenza.kadenza_tpf_main()
    except SystemExit:  # System exit is expected (required arguments missing)
        pass
    try:
        kadenza.kadenza_ffi_main()
    except SystemExit:
        pass

