A program which calculates the time-domain complex EM fields of an ultrashort tightly-focused Laguerre-Gaussian vortex beam of arbitrary orbital and radial modes in the perturbative regime.  This code is designed to calculate said fields to arbitrary perturbative order, as developed and generalizaed in our two related papers:

1) Development of the perturbative method

	Perturbative representation of ultrashort nonparaxial elegant Laguerre-Gaussian fields
	(https://journals.aps.org/pra/pdf/10.1103/PhysRevA.98.043820)

2) Generalization to arbitrary perturbative order

	Generalized Description of a Perturbative Nonparaxial Elegant Laguerre-Gaussian Phasor for Ultrashort Pulses in the Time Domain
	(https://journals.aps.org/pra/pdf/10.1103/PhysRevA.99.053832)


Configuration:
- Upgrade pip and supporting install tools: `pip install -U pip setuptools wheel`
- Install project requirements: `pip install -r requirements.txt`
- In constants.py, you can set various beam parameters such as the polarization, LG mode, spot size, and inital phase

Execution:
- Generate visualization in `pert_fields.ipynb`
- Generate exact fields at a particular space-time coordinate through `fieldsModule.makeFields`
