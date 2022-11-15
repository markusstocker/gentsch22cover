library(orkg)

orkg <- ORKG(host='https://sandbox.orkg.org/')

df <- orkg$resources$by_id('R205097')$as_dataframe()

print(df)