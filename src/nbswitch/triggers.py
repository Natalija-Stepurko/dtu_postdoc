"""What a pH trigger would actually have to work with.

Switching on pH means switching on protonation, and for a protein that means
histidine: it is the only side chain whose pKa sits near physiological pH, so
it is what every pH-responsive binder is built from. Whether that is usable
here is a question with an arithmetic answer, and the arithmetic is worth doing
before the idea is accepted or dismissed.

The Henderson-Hasselbalch equation gives the protonated fraction of a group:

    f = 1 / (1 + 10^(pH - pKa))

So the design question is not whether the pH differs between blood and a bite
site, but how much of the switching residue changes state across that
difference, and whether what remains is enough to act on.

**On the pKa.** A free histidine imidazole is near 6.0. In a protein it ranges
roughly 5.5 to 7.0 depending on environment, and engineering it is exactly what
a pH-switch design would do, so the single value used here is a reference point
and the sensitivity to it is reported alongside.
"""
from __future__ import annotations

#: Free histidine imidazole. See the module note: in a protein this moves.
HISTIDINE_PKA = 6.0

#: Arterial blood, the systemic state the construct has to stay masked in.
BLOOD_PH = (7.35, 7.45)

#: The range the case brief gives for a bite site. It is quoted from the brief
#: and is not a measurement made in envenomed tissue; see the page.
BITE_SITE_PH = (6.0, 6.5)


def fraction_protonated(ph: float, pka: float = HISTIDINE_PKA) -> float:
    """Henderson-Hasselbalch: the fraction of the group carrying a proton."""
    return 1.0 / (1.0 + 10.0 ** (ph - pka))


def switch_window(blood_ph=BLOOD_PH, site_ph=BITE_SITE_PH,
                  pka: float = HISTIDINE_PKA) -> dict:
    """How much of a histidine changes state between blood and a bite site.

    The conservative comparison is taken: the blood end that protonates most
    against the site end that protonates least, since a switch has to work at
    the worst combination rather than the most favourable one.
    """
    blood = {p: fraction_protonated(p, pka) for p in blood_ph}
    site = {p: fraction_protonated(p, pka) for p in site_ph}
    blood_worst = max(blood.values())        # most protonated while systemic
    site_worst = min(site.values())          # least protonated at the target
    return {
        "pka": pka,
        "blood_ph": list(blood_ph),
        "site_ph": list(site_ph),
        "blood_fraction": {str(k): round(v, 4) for k, v in blood.items()},
        "site_fraction": {str(k): round(v, 4) for k, v in site.items()},
        "blood_worst": round(blood_worst, 4),
        "site_worst": round(site_worst, 4),
        "fold_change": round(site_worst / blood_worst, 2),
        "unprotonated_at_site": round(1 - site_worst, 4),
    }


def sensitivity(pkas=(5.5, 6.0, 6.5, 7.0), blood_ph=BLOOD_PH,
                site_ph=BITE_SITE_PH) -> list[dict]:
    """The same comparison across the range a histidine pKa can take."""
    return [switch_window(blood_ph, site_ph, pka) for pka in pkas]
