# round((10.^(3:0.25:8)).^(1/3))
oilwater_hookean_DIMS := 10    12    15    18    22    26    32    38    46    56    68    83   100   121 #  147 178 215 #   261   316   383   464

$(foreach s,$(oilwater_hookean_DIMS),$(eval\
	$(call generic_osprey_helper,oilwater_hookean,oilwater_hookean-$sx$sx$s,$s,,True,) \
))
