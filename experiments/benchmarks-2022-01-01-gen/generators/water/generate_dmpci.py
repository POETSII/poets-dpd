import sys

sys.stderr.write(f"argv={sys.argv}\n")

dst=sys.argv[1]
scale=sys.argv[2]
steps=1
dt=0.01

frac_digits=4  # Positions have 4 fractional decimal digits of accuracy
mant_digits=5  # Vel and force have 5 relative decimal digits of accuracy

print(f"""dpd

Title	" Water "
Date    01/01/2
Comment	" Basic water analysis "

State	random

Bead  W1
      0.5
      25
      4.5

Bead  W2
      0.5
      25   25
      4.5  4.5

Polymer	Water1    0.99   " (W1) "
Polymer	Water2    0.01   " (W2) "

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

Command RoundBeadProperties 1 10 {frac_digits} {mant_digits}
Command ExportToPDPDWorldState 1 {dst}
Command StopNoSave 1
""")

