import sys

sys.stderr.write(f"argv={sys.argv}\n")

dst=sys.argv[1]
scale=sys.argv[2]
steps=1
dt=0.01

frac_digits=4  # Positions have 4 fractional decimal digits of accuracy
mant_digits=5  # Vel and force have 5 relative decimal digits of accuracy

print(f"""dpd

Title	" Water and oil phase separation "
Date    10/02/21
Comment	" 75:25 mixture of water and oil to see them phase separate. The cross interaction parameter
          is large to drive them to separate. dt10: I upped hookean bond strength from 65 to 100, as otherwise bonds were
          snapping early on and taking too long to smooth out.

          Note. If you edit the title above or this comment there must be at least one space between the quotes and the text. Blank lines are allowed.   "

State	random

Bead  W
      0.5
      15
      4.5

Bead  O
      0.5
      55    15
      4.5   4.5
      
<<<<<<< HEAD:experiments/benchmarks-2022-01-01/generators/oilwater_angles/generate_dmpci.py
Bond  O O  100.0  0.5
=======
Bond  O O  200.0  0.5
>>>>>>> c716c30c88d424543a372877c5110ec5463406d5:experiments/benchmarks-2022-01-01-gen/generators/oilwater_angles/generate_dmpci.py

BondPair O O O 5.0  0.0

Polymer	Water    0.75   " (W) "
Polymer Oil      0.25   " (O O O O O O) "

Box       {scale} {scale} {scale}       1  1  1
Density     3
Temp        1
RNGSeed     -7706
Lambda		0.5
Step		0.005
Time		  {steps*10}
SamplePeriod      {steps*10}
AnalysisPeriod    {steps*10}
DensityPeriod     {steps*10}
DisplayPeriod     {steps*10}
RestartPeriod     {steps*10}
Grid		1  1  1

Command SetTimeStepSize     {steps}     {dt}

Command OptimisePolymerOrderingForPDPD 1
Command RoundBeadProperties 1 10 {frac_digits} {mant_digits}
Command ExportToPDPDWorldState 1 {dst}
Command StopNoSave 1
""")

