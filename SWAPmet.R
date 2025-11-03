# -------------------------------
# نصب پکیج‌ها (در صورت نیاز)
# -------------------------------
# install.packages(c("nasapower", "dplyr", "lubridate"))

library(nasapower)
library(dplyr)
library(lubridate)

# -------------------------------
# 1. پارامترهای ایستگاه
# -------------------------------
lat <- 36.21   # عرض جغرافیایی
lon <- 57.68   # طول جغرافیایی
start_date <- "2023-04-21"
end_date   <- "2023-04-30"

# -------------------------------
# 2. دانلود داده از NASA POWER
# -------------------------------
vars <- c("T2M_MAX", "T2M_MIN", "PRECTOTCORR", "ALLSKY_SFC_SW_DWN", "RH2M", "WS2M")

np <- get_power(
  community = "AG",
  lonlat = c(lon, lat),
  pars = vars,
  temporal_api = "daily",
  dates = c(start_date, end_date)
)

df <- np %>%
  select(YEAR, DOY, YYYYMMDD, T2M_MAX, T2M_MIN, PRECTOTCORR, ALLSKY_SFC_SW_DWN, RH2M, WS2M) %>%
  rename(
    Year = YEAR,
    DayOfYear = DOY,
    Date = YYYYMMDD,
    Tmax = T2M_MAX,
    Tmin = T2M_MIN,
    Rain = PRECTOTCORR,
    Radn = ALLSKY_SFC_SW_DWN,
    RH = RH2M,
    Wind = WS2M
  ) %>%
  mutate(Date = as.Date(Date))

# -------------------------------
# 3. تابع محاسبه FAO-56 Penman–Monteith ET0
# -------------------------------
et0_fao56 <- function(Tmax, Tmin, Rs, u2, RH, lat, elev, doy){
  Tmean <- (Tmax + Tmin) / 2
  P <- 101.3 * ((293 - 0.0065 * elev) / 293)^5.26
  gamma <- 0.000665 * P
  delta <- 4098 * (0.6108 * exp((17.27 * Tmean) / (Tmean + 237.3))) / (Tmean + 237.3)^2
  es <- (0.6108 * exp((17.27 * Tmax) / (Tmax + 237.3)) + 0.6108 * exp((17.27 * Tmin) / (Tmin + 237.3))) / 2
  ea <- es * RH / 100
  
  Gsc <- 0.0820
  dr <- 1 + 0.033 * cos(2 * pi / 365 * doy)
  delta_s <- 0.409 * sin(2 * pi / 365 * doy - 1.39)
  phi <- lat * pi / 180
  ws <- acos(-tan(phi) * tan(delta_s))
  Ra <- (24 * 60 / pi) * Gsc * dr * (ws * sin(phi) * sin(delta_s) + cos(phi) * cos(delta_s) * sin(ws))
  Rso <- (0.75 + 2e-5 * elev) * Ra
  Rns <- (1 - 0.23) * Rs
  sigma <- 4.903e-9
  Rnl <- sigma * ((Tmax + 273.16)^4 + (Tmin + 273.16)^4)/2 *
    (0.34 - 0.14 * sqrt(ea)) * (1.35 * (Rs / Rso) - 0.35)
  Rn <- Rns - Rnl
  G <- 0
  
  ET0 <- (0.408 * delta * (Rn - G) + gamma * (900 / (Tmean + 273)) * u2 * (es - ea)) /
    (delta + gamma * (1 + 0.34 * u2))
  return(ET0)
}

# -------------------------------
# 4. محاسبه ET0 برای هر روز
# -------------------------------
elev <- 985  # ارتفاع از سطح دریا (متر)
df <- df %>%
  mutate(doy = yday(Date)) %>%
  rowwise() %>%
  mutate(ETref = et0_fao56(Tmax, Tmin, Radn, Wind, RH, lat, elev, doy)) %>%
  ungroup()

# -------------------------------
# 5. تبدیل به فرمت SWAP
# -------------------------------
swap_df <- df %>%
  mutate(
    Station = "'Sabzevar'",
    MM = month(Date),
    DD = day(Date),
    YYYY = year(Date)
  ) %>%
  select(Station, DD, MM, YYYY, Radn, Tmin, Tmax, RH, Wind, Rain, ETref)

# -------------------------------
# 6. ذخیره فایل خروجی
# -------------------------------
write.table(swap_df, "SWAPmet_Sabzevar.txt", sep = "\t", row.names = FALSE, quote = FALSE)

cat("✅ فایل SWAPmet_Sabzevar.txt با موفقیت ایجاد شد!\n")
