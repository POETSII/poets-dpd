import sys

sys.stderr.write(f"argv={sys.argv}\n")

dst=sys.argv[1]
scale=sys.argv[2]
steps=1
dt=0.02

frac_digits=4  # Positions have 4 fractional decimal digits of accuracy
mant_digits=5  # Vel and force have 5 relative decimal digits of accuracy

print(f"""dpd

Title	" Water and oil phase separation "
Date    10/02/21
Comment	" 75:25 mixture of water and oil to see them phase separate. The cross interaction parameter
          is large to drive them to separate.

          Note. If you edit the title above or this comment there must be at least one space between the quotes and the text. Blank lines are allowed.   "

State	random

Bead  W
      0.5
      25
      4.5

Bead  O
      0.5
      75    25
      4.5   4.5

Polymer	Water    0.75   " (W) "
Polymer Oil      0.25   " (O) "

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

