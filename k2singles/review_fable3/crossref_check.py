"""Fable review 3: verify bibliography entries against Crossref (journal, volume, page/article, first author, year)."""
import json, urllib.request, urllib.parse
items = {
 "Armstrong2016": ("10.1093/mnras/stv2836", "MNRAS 456 2260"),
 "Caldwell2020": ("10.3847/2515-5172/abc9b3", "RNAAS 4 201"),
 "Arenou2023": ("10.1051/0004-6361/202243782", "A&A 674 A34"),
 "Gaia2023": ("10.1051/0004-6361/202243940", "A&A 674 A1"),
 "Ginsburg2019": ("10.3847/1538-3881/aafc33", "AJ 157 98"),
 "Guo2025": ("10.3847/1538-4365/adc38a", "ApJS 278 46"),
 "Howell2014": ("10.1086/676406", "PASP 126 398"),
 "Huang2020": ("10.3847/2515-5172/abca2e", "RNAAS 4 204"),
 "IJspeert2024": ("10.1051/0004-6361/202450505", "A&A 691 A242"),
 "Kovacs2020": ("10.1051/0004-6361/202038725", "A&A 643 A169"),
 "Kovalev2024": ("10.1093/mnras/stad3105", "MNRAS 527 521"),
 "LaCourse2018": ("10.3847/2515-5172/aaad61", "RNAAS 2 28"),
 "Libralato2016": ("10.1093/mnras/stw1932", "MNRAS 463 1780"),
 "Sandquist2021": ("10.3847/1538-3881/abca8d", "AJ 161 59"),
 "Valizadegan2025": ("10.3847/1538-3881/ae0a5c", "AJ 170 287"),
}
def q(url):
    req = urllib.request.Request(url, headers={"User-Agent": "duo-transits review (mailto:jon.mcbee@gmail.com)"})
    return json.load(urllib.request.urlopen(req, timeout=60))
for key, (doi, expect) in items.items():
    try:
        m = q(f"https://api.crossref.org/works/{urllib.parse.quote(doi)}")["message"]
        au = m.get("author", [])
        a0 = (au[0].get("family", "") + " " + au[0].get("given", "")) if au else "?"
        yr = (m.get("published-print") or m.get("published-online") or m.get("issued"))["date-parts"][0][0]
        print(f"{key:16s} DOI ok | {m.get('container-title',['?'])[0][:40]:40s} vol {m.get('volume')} page {m.get('page') or m.get('article-number')} | {yr} | {a0} (+{len(au)-1}) | {m.get('title',[''])[0][:90]} | expect {expect}")
    except Exception as ex:
        # fall back: bibliographic query
        try:
            title = {"Guo2025": "Guo LAMOST SB2 orbit ApJS 2025", "Valizadegan2025": "Valizadegan ExoMiner 2025 AJ", "IJspeert2024": "IJspeert eclipsing binaries TESS OBAF 2024", "Kovacs2020": "Kovacs 2020 A&A single transit K2", "Kovalev2024": "Kovalev LAMOST MRS SB2 orbits 2024 MNRAS"}.get(key, key)
            r = q("https://api.crossref.org/works?rows=3&query.bibliographic=" + urllib.parse.quote(title))["message"]["items"]
            print(f"{key:16s} DOI {doi} FAILED ({ex}); bibliographic search:")
            for m in r:
                au = m.get("author", []); a0 = au[0].get("family", "") if au else "?"
                print(f"      {m.get('DOI')} | {m.get('container-title',['?'])[0][:30]} vol {m.get('volume')} page {m.get('page') or m.get('article-number')} | {a0} | {m.get('title',[''])[0][:80]}")
        except Exception as ex2:
            print(f"{key:16s} FAILED both: {ex} / {ex2}")
