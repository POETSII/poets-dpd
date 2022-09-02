import sys

sys.stderr.write(f"argv={sys.argv}\n")

dst=sys.argv[1]
scale=sys.argv[2]
steps=1
dt=0.02

frac_digits=4  # Positions have 4 fractional decimal digits of accuracy
mant_digits=5  # Vel and force have 5 relative decimal digits of accuracy

print(f"""dpd

Title	" Water/Amphiphile bilayer "
Date    10/02/21
Comment	" The conservative repulsion parameters are taken from Grafmueller et al. Biophysical Journal 96:2658 (2009).
          Warning: this run will take a long time if the default length is left at 200,000 time steps. But that is
          the only way to get statistical accuracy for the surface tension.   "

State	lamella
       	Polymer			Lipid
        Normal			0 0 1
        Centre			0.5
        Thickness		5.0
        Linearise		1
       	UpperFraction	0.5
       	Polymerise		0

Bead  H
      0.5
      30	
      4.5	

Bead  T
      0.5
      35  10
      4.5 4.5

Bead  W
      0.5
      30	75	25
      4.5	4.5	4.5


Bond    H H  128	0.5
Bond    H T  128	0.5
Bond    T T  128	0.5

BondPair	H T T	15.0	0.0
BondPair	T T T	15.0	0.0

Polymer	Water  0.974223   " (W) "
Polymer	Lipid  0.025777   " (H H (* (T T T T T T)) H T T T T T T) "

Box       {scale} {scale} 32       1  1  1
Density     3
Temp        1
RNGSeed     -7706
Lambda		0.5
Step		0.005
Time		      {10000}
SamplePeriod      {10}
AnalysisPeriod    {10}
DensityPeriod     {10}
DisplayPeriod     {10}
RestartPeriod     {10}
Grid		1  1  1

Command OptimisePolymerOrderingForPDPD 1
Command RoundBeadProperties 1 10 {frac_digits} {mant_digits}
Command ExportToPDPDWorldState 1 {dst}

Command SetCurrentStateDefaultFormat 1 Paraview

Command	ToggleBeadDisplay       1    W
Command SetCurrentStateCamera   1    0.5 -0.5 -0.5  0.5 0.5 0.5

Command StopNoSave 1
""")



