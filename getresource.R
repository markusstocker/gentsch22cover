library(orkg)

orkg <- ORKG(host='https://sandbox.orkg.org/')

df <- orkg$resources$by_id('R202839')$as_dataframe()

print(df)