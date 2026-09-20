
def check(a, b):
    assert abs(a - b)/abs(b) < 1e-6, f"{a} != {b}"


def check_system(u):

    # check(u.oz, u.lbm/16)

    check(u.lbf/u.slug,   u.ft / u.sec**2)
    check(u.pdl,   0.138254954376*u.N)
    check(u.slug, 32.17405*u.lbm)
    check(u.slug,  14.59390*u.kg)
    check(u.lbm, 0.45359237*u.kg)  # International avoirdupois pound


    check(u.lbf/u.slug,   u.ft / u.sec**2)
    check(u.pdl,   0.138254954376*u.N)
    check(u.slug, 32.17405*u.lbm)
    check(u.slug,  14.59390*u.kg)
    check(u.lbm, 0.45359237*u.kg)  # International avoirdupois pound

    check(u.minute,   60*u.sec)
    check(u.hr,     3600*u.sec)
    check(u.day,      24*u.hr)
    check(u.yr,   365.25*u.day)

    check(u.inch,   25.4*u.mm)     # International inch
    check(u.ft,       12*u.inch)
    check(u.yd,        3*u.ft)
    check(u.mi,     5280*u.ft)
    check(u.mi, 1.609344*u.km)
    check(u.km,     1000*u.m)
    check(u.m,       100*u.cm)

    check(u.N,   u.kg*u.m/u.sec**2)
    check(u.dyn, u.gm*u.cm/u.sec**2)
    check(u.pdl, u.lbm*u.ft/u.sec**2)
    check(u.kN,  1000*u.N)
    check(u.dyn,  1e-5*u.N)
    check(u.lbf, 4.4482216152605*u.N)
    check(u.kip, 1000*u.lbf)
    check(u.kgf, 9.80665*u.N)

    check(u.gm,     1e-3*u.kg)
    check(u.tonne,  1000*u.kg)
    check(u.slinch, u.lbf*u.sec**2/u.inch)
    check(u.slinch,   12*u.slug)
    check(u.slinch, 175.1268*u.kg)

    # Standard gravity
    check(u.gravity,  9.80665*u.m/u.sec**2)
    check(u.gravity, 32.17405*u.ft/u.sec**2)
    check(u.gravity, 386.0886*u.inch/u.sec**2)
    check(u.lbf, u.lbm*u.gravity)
    check(u.kgf, u.kg*u.gravity)

    # check(180*u.deg, u.pi*u.rad)


def test_units():
    from xara.units import si, iks, ips, fks, fps

    check_system(si)
    check_system(iks)
    check_system(ips)
    check_system(fks)
    check_system(fps)

