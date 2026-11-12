# Criterio de vecindad para SESCC-pairs (fijado el 11-11-2026, antes de aplicarlo al Almagesto)

- Residuos de longitud en unidades angulares (`--cosb`): los errores de longitud de
  los tres catálogos son angulares (dispersión de Δλ·cosβ constante con β; check_noise.py).
- Vecindad por separación angular: `--sep 15`, centro del intervalo 10–30° en el que
  Tycho y Ulugh Beg recuperan su época (validacion_cosb_sep.log).
- Normalización de pares sin cambios: `--pairnorm 0.5`.
- Se publicará la sensibilidad completa (10°, 15°, 20°, 30° y la caja Δλ<30° anterior).

## Nota sobre el registro
Este fichero se escribió antes de aplicar el criterio al Almagesto (22:24, 11-11-2026),
pero por un error en la orden `git add` no entró en el commit d774e02 y se confirmó
después, cuando ya se habían visto los resultados del Almagesto con `--sep 10, 15, 20
y 30 --cosb` (almagest_cosb_sep.log: -134, -274, -223, -223). Por tanto, este commit
no demuestra un registro previo. El criterio se justifica por la estructura del ruido
(check_noise.py) y por los catálogos de validación (validacion_cosb_sep.log), y se
publica la sensibilidad completa.

## Corrección posterior (12-11-2026)
Los sintéticos de longitud con `--sep` mostraban un atractor espurio en +232. Causa: pares
que rodean el polo de la eclíptica (|Δλ| cerca de 180°), cuyo residuo saltaba 360° y para
los que Δλ·cosβ no aproxima la separación. Corrección: residuo reducido a ±180° y exclusión
de los pares con |Δλ|·cos(β medio) mayor que la separación máxima (37 de 8492 pares en el
Almagesto a 15°). Con el código corregido, el criterio (`--sep 15 --cosb`) sigue
recuperando las épocas de Tycho y Ulugh Beg (sensibilidad_v2c.log), y las estimaciones del
Almagesto apenas cambian. La caja |Δλ|<30° anterior da los mismos resultados que antes.

## Signed pair weights (added after the v2c runs)
Eq. (2) defines the pair weight with signed angular speeds,
|mu_i cos(b_i) - mu_j cos(b_j)|. Up to v2c, `pmotion_lon` returned absolute
values, so the weight was ||v_i| - |v_j||. The option `--signed` implements
Eq. (2); the default keeps the old behaviour for reproducibility.
The criterion (`--sep 15 --cosb`) was NOT re-chosen: with `--signed` it gives
Tycho 1572 and Ulugh Beg 1457, as before (sensibilidad_v2d.log).
Almagest: -274 -> -240.
