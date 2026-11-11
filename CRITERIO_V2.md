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
