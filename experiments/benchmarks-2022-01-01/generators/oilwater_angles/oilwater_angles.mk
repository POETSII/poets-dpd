# round((10.^(3:0.25:8)).^(1/3))
OILWATER_ANGLES_DIMS := 10    12    15    18    22    26    32    38    46    56    68    83   100   121  # 147 178 215 #   261   316   383   464

$(foreach s,$(OILWATER_ANGLES_DIMS),$(eval\
	$(call generic_osprey_helper,oilwater_angles,oilwater_angles-$sx$sx$s,$s,,True,True) \
))
