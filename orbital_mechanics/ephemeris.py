from skyfield.api import load
from skyfield.elementslib import osculating_elements_of

ts = load.timescale()
t = ts.utc(2018, 4, 22, range(0, 25))

planets = load('de421.bsp')
earth = planets['earth']
moon = planets['moon']

position = (moon - earth).at(t)
elements = osculating_elements_of(position)

a = '321'
