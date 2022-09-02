# round((10.^(3:0.25:8)).^(1/3))
oilwater_hookean_DIMS := 10    12    15    18    22    26    32    38    46    56    68    83   100   #121  147 178 215 #   261   316   383   464

$(foreach s,$(oilwater_hookean_DIMS),$(eval\
<<<<<<< HEAD:experiments/benchmarks-2022-01-01/generators/oilwater_hookean/oilwater_hookean.mk
	$(call generic_osprey_helper,oilwater_hookean,oilwater_hookean-$sx$sx$s,$s,2000,True,) \
=======
	$(call generic_osprey_helper,oilwater_hookean,oilwater_hookean-$sx$sx$s,$s,1000,True,) \
>>>>>>> c716c30c88d424543a372877c5110ec5463406d5:experiments/benchmarks-2022-01-01-gen/generators/oilwater_hookean/oilwater_hookean.mk
))
