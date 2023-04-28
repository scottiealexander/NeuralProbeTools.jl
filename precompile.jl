using OEphysCSD
basedir = "/home/scottie/data/oephys/recording3"
h, _, _ = OEphysCSD.run(basedir, bad_channels=OEphysCSD.REC3_BAD)
close(h)
h, _, _ = OEphysCSD.run(basedir)
close(h)
