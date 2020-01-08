# `cdips_followup`

Given planet candidates, turn them into planets.

Related development also happens in the `/cdips/` repo.

This repo is for:

* analyses of individual objects of interest.
* ephemeris tracking, via `ephemerides.csv`
* candidate tracking, via `candidates.csv`
* LCOGT followup requests

__install__

`conda env create -f environment.yml -n cdips_followup`

`python setup.py develop`

__update__

`conda env export > environment.yml`
