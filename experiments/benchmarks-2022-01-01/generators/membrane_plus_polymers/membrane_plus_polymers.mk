# round((10.^(3:0.25:8)/48).^(1/2))
membrane_plus_polymers_DIMS := 5      6      8     11     14     19     26     34     46     61     81    108    144    192    257    342    456  #  609  812   1082   1443

$(foreach s,$(membrane_plus_polymers_DIMS),$(eval \
	$(call generic_osprey_helper,membrane_plus_polymers,membrane_plus_polymers-$sx$sx48,$s,1000,True,True)\
))


