#############
ParticleLoad
#############

**Description:** Particle property which carries a load or mass, which can be used to calculate spatial and temporal statistics of this load or mass

**class_name:** oceantracker.particle_properties.load_carrying.ParticleLoad

**File:** oceantracker/particle_properties/load_carrying.py

**Inheritance:** _BasePropertyInfo> ParticleProperty> ParticleLoad


Parameters:
************

	* ``class_name`` :   ``<class 'str'>``   *<optional>*
		Description: Class name as string A.B.C, used to import this class from python path

		- default: ``None``

	* ``description`` :   ``<class 'str'>``   *<optional>*
		- default: ``None``

	* ``dtype`` :   ``<class 'numpy.dtype'>``   *<optional>*
		- default: ``<class 'numpy.float64'>``

	* ``fill_value`` :   ``[<class 'int'>, <class 'float'>]``   *<optional>*
		- default: ``None``

	* ``initial_value`` :   ``<class 'float'>``   *<optional>*
		Description: Value of the particle property when it is released

		- default: ``1.0``

	* ``prop_dim3`` :   ``<class 'int'>``   *<optional>*
		- default: ``1``
		- min: ``1``

	* ``release_group_parameters``: nested parameter dictionary
	* ``time_varying`` :   ``<class 'bool'>``   *<optional>*
		- default: ``True``
		- possible_values: ``[True, False]``

	* ``type`` :   ``<class 'str'>``   *<optional>*
		Description: type of particle property, used to manage how to update particle property

		- default: ``user``
		- possible_values: ``['manual_update', 'from_fields', 'user']``

	* ``update`` :   ``<class 'bool'>``   *<optional>*
		- default: ``True``
		- possible_values: ``[True, False]``

	* ``user_note`` :   ``<class 'str'>``   *<optional>*
		- default: ``None``

	* ``vector_dim`` :   ``<class 'int'>``   *<optional>*
		- default: ``1``
		- min: ``1``

	* ``write`` :   ``<class 'bool'>``   *<optional>*
		Description: Write particle property to tracks or event files file

		- default: ``True``
		- possible_values: ``[True, False]``

