library(orkg)

orkg <- ORKG(host='https://sandbox.orkg.org/')

r <- orkg$resources$by_id('R205097')

print(r)