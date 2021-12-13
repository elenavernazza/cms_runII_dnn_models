# cms_runII_dnn_models

Repo for models used in CMS Run II hh->bbtautau search

For use with:

- Computation interface: [cms_hh_proc_interface](https://github.com/GilesStrong/cms_hh_proc_interface)
- Model interface: [cms_hh_tf_inference](https://github.com/GilesStrong/cms_hh_tf_inference)

# Non-res GluGlu

|Model Name|Jet Ordering|Mass Cut|Interface Version|
|---|---|---|---|
|`2020-06-29-0`|DF|Default|V3.0 or V4.0|
|`2020-06-29-1`|DF|Setting 1|V4.0|
|`2020-06-30-0`|HHBTag|Default|V3.0 or V4.0|
|`2020-06-30-1`|HHBTag|Setting 2|V4.0|
|`2020-07-31-0`|HHBTag|Setting 2|V3.0 or V4.0|

# ARC checks

## zz_bbtt

|Model Name|Jet Ordering|Mass Cut|Interface Version|
|---|---|---|---|
|`2021-11-22-0`|HHBTag|Setting ZZ|V3.0 or V4.0|

## zh_bbtt

|Model Name|Jet Ordering|Mass Cut|Interface Version|
|---|---|---|---|
|`2021-11-22-0`|HHBTag|Setting ZH|V3.0 or V4.0|

# Res GluGlu

|Model Name|Jet Ordering|Mass Cut|Parameterised?|Interface Version|
|---|---|---|---|---|
|`2020-08-01-0`|HHBTag|Setting 2|Yes|V3.0 or V4.0|
|`2020-08-03-0`|HHBTag|Setting 2|No|V3.0 or V4.0|

# Mass Cuts

### Non-boosted

|Mass cut|m_tt offset|m_tt res.|m_bb offset|m_bb res.|
|---|---|---|---|---|
|Default|116|35|111|45|
|Setting 1|127|50|188|193|
|Setting 2|129|53|169|145|
|Setting ZZ|115|182|105|280|
|Setting ZH|140|147|80|228|

### Boosted window

|Mass cut|m_tt min.|m_tt max|m_bb min.|m_bb max.|
|---|---|---|---|---|
|Default|80|152|90|160|

### Boosted ellipse

|Mass cut|m_tt offset|m_tt res.|m_bb offset|m_bb res.|
|---|---|---|---|---|
|Setting 1|126|81|181|110|
|Setting 2|128|60|159|94|
|Setting ZZ|97|96|48|33|
|Setting ZH|112|181|84|115|