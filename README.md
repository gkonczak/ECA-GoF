---
title: "Untitled"
author: "Grzegorz Kończak"
date: "2026-02-06"
output: html_document
---

# eca_gof.test: Empty Cell Area Goodness-of-Fit Test

## Overview

**eca_gof.test** is an **R** function implementing the **Empty Cell Area (ECA) test** — a nonparametric goodness-of-fit test that evaluates whether a sample comes from a specified theoretical distribution. The test works by measuring the proportion of the probability space that remains uncovered by neighborhoods around sample points.

## Usage

eca.test(x, dist = "norm", params = list(mean = 0, sd = 1),
         delta = 1 / (2 * length(x)), alpha = 0.05, n.sim = 1000)
         
         
## Arguments


x	numeric vector	—	Sample observations

dist	character	"norm"	Distribution family ("norm", "exp", "unif", "gamma", "lnorm", …)

params	named list	list(mean=0, sd=1)	Parameters passed to d/p/q/r functions

delta	numeric	1/(2*n)	Half-width of probability neighborhood

alpha	numeric	0.05	Significance level

n.sim	integer	1000	Number of Monte Carlo simulations

## Values

statistic	Observed ECA test statistic

p.value	Monte Carlo p-value

alpha	Significance level used

decision	"REJECT H0" or "DO NOT REJECT H0"

method	Description of the test

data.name	Name of the input data

