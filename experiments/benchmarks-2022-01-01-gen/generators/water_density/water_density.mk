WATER_DENSITY_DIMS := 25 35 45 55
WATER_DENSITY_DENSITIES := 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 

$(foreach s,$(WATER_DENSITY_DIMS),$(foreach d,$(WATER_DENSITY_DENSITIES),$(eval\
	$(call generic_osprey_helper,water_density,water_density-$sx$sx$s-d$d,$s $d) \
)))
