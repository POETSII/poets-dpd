# round((10.^(3:0.25:8)/32).^(1/2))
bilayer_DIMS := 6      7     10     13     18     24     31     42     56     75     99    133    177    236    314    419    559 #   745    994   1326   1768

$(foreach s,$(bilayer_DIMS),$(eval \
	$(call generic_osprey_helper,bilayer,bilayer-$sx$sx32,$s,1000,True,True)\
))


