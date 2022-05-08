import sys

sys.stderr.write(f"argv={sys.argv}\n")

dst=sys.argv[1]
scale=sys.argv[2]

frac_digits=4  # Positions have 4 fractional decimal digits of accuracy
mant_digits=5  # Vel and force have 5 relative decimal digits of accuracy

print(f"""dpd

Title	" LLPS at a membrane "
Date    08/10/21
Comment  " Two types of H3(T4)2 lipid with the polymers attracted to the head groups of type B only.     "

State	compositelamella
        Polymers		Lipid CoLipid
        Normal			0 0 1
        Centre			0.75
        Thickness		5.0
        Linearise		1
        UpperFraction	0.6  0.0
        Patches         0  0
        Polymerise		0

Bead  W
      0.5
      25
      4.5

Bead  E
      0.5
      25    5
      4.5   4.5

Bead  B
      0.5
      23    25    25
      4.5   4.5   4.5

Bead  HA
      0.5
      30    30   30   30
      4.5   4.5  4.5  4.5

Bead  TA
      0.5
      75    35   35   35   10
      4.5   4.5  4.5  4.5  4.5

Bead  HB
      0.5
      30    5    30   30   35   30
      4.5   4.5  4.5  4.5  4.5  45

Bead  TB
      0.5
      75   35   35   35   10   35   10
      4.5  4.5  4.5  4.5  4.5  4.5  4.5



Bond	E  E    128  0.5
Bond	B  B    128  0.5
Bond	E  B    128  0.5
Bond    HA HA   128  0.5
Bond    HA TA   128  0.5
Bond    TA TA   128  0.5

Bond    HB HB   128  0.5
Bond    HB TB   128  0.5
Bond    TB TB   128  0.5


BondPair  B B B      5.0   0.0
BondPair  HA TA TA  15.0   0.0
BondPair  TA TA TA  15.0   0.0
BondPair  HB TB TB  15.0   0.0
BondPair  TB TB TB  15.0   0.0


Polymer Water    0.984    " (W) "
Polymer Lipid    0.0118   " (HA HA HA (* (TA TA TA TA)) HA TA TA TA TA) "
Polymer CoLipid  0.0002   " (HB HB HB (* (TB TB TB TB)) HB TB TB TB TB) "
Polymer Rod      0.004    " (E E (* E) (* E) (8 B) E (* E) (* E) E)  "

Box       {scale} {scale}  48       1  1  1
Density     3
Temp        1
RNGSeed     -7706
Lambda		0.5
Step		0.01
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



