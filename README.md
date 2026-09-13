# DustVel

Chang zarrachalarining atmosferadagi harakat tezliklarini va traektoriyasini hisoblovchi dasturiy vosita.

[![Python](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

## Tavsif

DustVel shamol eroziyasi va texnologik jarayonlar natijasida atmosferaga ko'tarilgan chang zarrachalarining harakat tezliklarini hisoblaydi. Mavjud ko'pchilik modellardan farqli ravishda, dastur zarracha tezligini shamol tezligiga tenglashtirmaydi, balki uning inersion xossalarini hisobga olgan holda alohida aniqlaydi.

Dastur Buxoro davlat universitetida bajarilgan dissertatsiya ishi doirasida yaratilgan.

## Matematik asos

Zarracha harakati Nyutonning ikkinchi qonuni va Stoks qarshilik kuchi asosida tavsiflanadi:

```
dv_x/dt = beta (u_x - v_x)
dv_y/dt = beta (u_y - v_y)
dv_z/dt = -beta v_z - g
```

bu yerda `beta = 9 mu_a / (2 rho_p r^2)`, `tau = 1/beta` — relaksatsiya vaqti.

## Sonli algoritm

Tizim qattiq (stiff) bo'lgani uchun oshkor sxemalar qo'llanilmaydi. Makroqadam `dt` uchta teng qismga bo'linadi (`dtau = dt/3`) va har bir makroqadamda ketma-ket uchta bosqich bajariladi:

| Bosqich | Sxema | Aniqlik | Turg'unlik |
|---------|-------|---------|------------|
| `n -> n+1/3` | Krank–Nikolson | O(dtau²) | A-turg'un |
| `n+1/3 -> n+2/3` | BDF2 | O(dtau²) | A- va L-turg'un |
| `n+2/3 -> n+1` | BDF2 | O(dtau²) | A- va L-turg'un |

Bunday tuzilish algoritmni o'z-o'zidan boshlanuvchi (self-starting) qiladi. Krank–Nikolson sxemasi L-turg'un bo'lmagani uchun qattiq rejimda tebranishlar hosil qilishi mumkin, keyingi ikkita L-turg'un BDF2 bosqichi esa ularni so'ndiradi.

Traektoriya makroqadam ichidagi to'rtta tugun bo'yicha trapetsiya formulasi bilan integrallanadi, bu esa umumiy ikkinchi tartibli aniqlikni saqlaydi.

## O'rnatish

```bash
git clone https://github.com/<foydalanuvchi>/DustVel.git
cd DustVel
pip install -r requirements.txt
```

## Foydalanish

### Buyruqlar qatori

```bash
python -m dustvel --ux 7 --h0 25 --vaqt 200 --taqqoslash
```

Asosiy parametrlar:

| Parametr | Tavsif | Standart |
|----------|--------|----------|
| `--ux` | shamol tezligi, m/s | 5.0 |
| `--h0` | chiqarish balandligi, m | 20.0 |
| `--vaqt` | hisoblash davri, s | 200.0 |
| `--dt` | vaqt qadami, s | 1e-3 |
| `-d` | zarracha diametri, m | — |
| `--rho-p` | zarracha zichligi, kg/m³ | 2500 |
| `--harorat` | havo harorati, K | 293.15 |
| `--fraksiya-fayl` | fraksiyalar JSON fayli | — |
| `--chiqish` | natijalarni CSV ga saqlash | — |
| `--taqqoslash` | soddalashtirilgan model bilan taqqoslash | — |

### Grafik interfeys

```bash
python -m dustvel --gui
```

Interfeysda shamol tezligi, chiqarish balandligi, havo harorati va vaqt qadami surgichlar orqali o'zgartiriladi. Uchta ko'rinish rejimi mavjud: tezliklar, traektoriya va soddalashtirilgan model bilan taqqoslash.

### Kutubxona sifatida

```python
from dustvel import Fraksiya, Muhit, Hisoblagich

muhit = Muhit(ux=5.0, T_muhit=293.15)
h = Hisoblagich(muhit)
fr = Fraksiya("50 um", 50e-6, 2500.0)

print(f"tau = {h.tau(fr)*1e3:.3f} ms")
print(f"v_t = {h.chokish_tezligi(fr):.4f} m/s")

n = h.hisobla(fr, T=200.0, dt=1e-3, H0=20.0)
print(f"yerga yetish: t = {n.yer_t:.1f} s, x = {n.yer_x:.1f} m")
```

## Testlar

```bash
python tests/test_algoritm.py
```

Testlar quyidagilarni tekshiradi: relaksatsiya vaqti formulasi, cho'kish tezligi, yaqinlashish tartibi (O(dt²) ekani), katta vaqt qadamida turg'unlik va impulsli ko'tarilish sharti.

## Namunaviy natijalar

`ux = 5` m/s, `H0 = 20` m, `dt = 1e-3` s da:

| Fraksiya | d, µm | rho_p, kg/m³ | tau, ms | v_t, m/s | Re_p |
|----------|-------|--------------|---------|----------|------|
| PM2.5 | 2.5 | 1500 | 0.029 | 2.82e-04 | 0.000 |
| PM10 | 10 | 2000 | 0.613 | 6.01e-03 | 0.004 |
| Yirik | 50 | 2500 | 19.148 | 1.88e-01 | 0.635 |

50 µm li zarracha 20 m balandlikdan 106 s ichida yerga yetadi va 532 m masofani bosib o'tadi.

## Cheklovlar

- Qarshilik kuchi Stoks yaqinlashuvida hisoblanadi. `Re_p > 1` bo'lganda dastur ogohlantirish beradi; 60–80 µm dan yirik zarrachalar uchun nochiziqli qarshilikka tuzatma talab qilinadi.
- Submikron zarrachalar uchun Kanningem sirpanish tuzatmasi hisobga olinmagan.
- Turbulent diffuziya va konvektiv oqimlar qaralmaydi; ular ko'chish-diffuziya tenglamasi bilan bog'langanda hisobga olinadi.

## Havolalar

Ravshanov N., Shafiyev T.R., Bobojonova M.A., Kobilova D.O. Issledovanie klyuchevyx faktorov, vliyayushchix na rasprostranenie pylevyx chastits v atmosfere // Problemy vychislitel'noy i prikladnoy matematiki. — 2026. — № 3(73). — S. 2–18.

## Litsenziya

MIT — [LICENSE](LICENSE) fayliga qarang.
