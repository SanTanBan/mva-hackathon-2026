#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Narrated slide video for the Track 2 pitch, MVA Hackathon 2026.

Kaggle script kernel `santanubanerjee9/mva-pitch-video` (private, CPU only, internet on).

Public content only. The spoken text and the slide wording are embedded below and come
from reports/track2/pitch_script.md, with slide wording checked against the public
reports/track2/track2_report.md. No patient data, no names, dates or places beyond "UK",
no variant coordinates, no doses and no treatment advice. Candidates are hypotheses for
expert follow-up, never recommendations.

Steps, all inside the kernel:
  1. install the TTS: Kokoro-82M via pip `kokoro` + `misaki[en]` and espeak-ng via apt.
     Kokoro runs in a worker process, in the kernel's Python if it installs there (on
     Python 3.13 only with --ignore-requires-python), else in a uv-made Python 3.12 venv,
     else the narration falls back to Piper,
  2. render 1920x1080 slides with PIL (drawn at 2x, then downsampled),
  3. narrate each script section (fallback: Piper TTS); slide timing follows the audio,
  4. assemble an H.264 + AAC MP4 with ffmpeg,
  5. write summary.json, captions.srt, narration.txt, phonemes.txt, 4 thumbnails and a
     contact sheet so the result can be checked without playing it.

Outputs (/kaggle/working): mva_track2_pitch.mp4, summary.json, captions.srt, narration.txt,
phonemes.txt, thumb_1..4_*.png, contact_sheet.png, slides/slide_*.png.
"""

import json
import math
import os
import platform
import re
import shutil
import string
import subprocess
import sys
import textwrap
import time
import traceback
import wave
from datetime import datetime, timezone
from pathlib import Path

T0 = time.time()

# The source must reach Kaggle as UTF-8 (IPA phonemes, arrows, quotes below).
assert 'ˈɛ→·' == 'ˈɛ→·', 'kernel source was not decoded as UTF-8'

os.environ.setdefault('HF_HUB_DISABLE_PROGRESS_BARS', '1')
os.environ.setdefault('TOKENIZERS_PARALLELISM', 'false')
os.environ.setdefault('PYTHONWARNINGS', 'ignore')

# ----------------------------------------------------------------------------- settings
OUT = Path('/kaggle/working') if Path('/kaggle/working').is_dir() else Path.cwd() / 'out'
WORK = Path('/tmp/mva_pitch')
VIDEO_NAME = 'mva_track2_pitch.mp4'
W, H, SS, FPS = 1920, 1080, 2, 30
REPO_URL = ''  # public repository URL for the closing slide; empty = no link is shown

VOICE, LANG = 'af_heart', 'a'          # Kokoro: American English, highest-rated voice
PIPER_VOICE = 'en_US-lessac-high'      # fallback only

# Timeline, in seconds. Slide durations follow the narration of each section.
TIMING = dict(title_hold=4.0, lead=0.5, tail=0.7, close_hold=12.0, sentence_pause=0.30)
TARGET_TOTAL, MIN_TOTAL, MAX_TOTAL, HARD_MAX = 170.0, 160.0, 175.0, 179.0
TARGET_LUFS = -16.0                    # integrated loudness of the finished video (EBU R128)
SPEED_MIN, SPEED_MAX = 0.92, 1.10      # never faster than 1.1x; words are never cut
CONTENT_LIMIT = 940                    # slide content must end well above the footer rule (y=986)
FONT_SCALES = [1.16, 1.12, 1.08, 1.04, 1.0, 0.96, 0.92, 0.88, 0.84, 0.8]  # largest that fits wins

FOOTER = 'Research hypotheses for expert follow-up. Not medical advice.'
KICKER = 'MVA HACKATHON 2026 · TRACK 2'

# ----------------------------------------------------------------------------- script
# Spoken lines of reports/track2/pitch_script.md, verbatim (markdown italics removed; the
# bracketed [v2: ...] production note on slide 4 is not spoken).
NARRATION = [
    "Every drug in this pitch has to survive one rule first. If it shortens mitosis in the "
    "patient's own cells, it is out, however good its other results look. It comes straight "
    "from the biology of this disease, and it rules out several obvious drug classes before "
    "we start.",

    "MVA type 1 is a dosage disease. One copy of BUB1B makes no protein, the other makes too "
    "little, and adding BUBR1 back rescues patient cells. So the goal is more BUBR1, not a "
    "stronger kinase: BUBR1 is a pseudokinase. Its level is controlled by an acetylation "
    "cycle, and mice that cannot acetylate it get shorter mitoses and look like MVA. Shorter "
    "mitosis is the fingerprint of the disease. That is the veto.",

    "We checked five trial registries, PubMed, ChEMBL and Open Targets. None lists a trial "
    "for MVA. The only positive tractability signals are traps: one comes from the kinase "
    "family this protein no longer behaves like, and the other says it can be degraded, "
    "which is the opposite of what a dosage disease needs.",

    "We rank five hypotheses by how cheaply the first experiment could disprove them. NAD+ "
    "precursors and HDAC2/3 inhibition both aim to protect BUBR1 itself. Metformin could sit "
    "inside the kidney surveillance already recommended in MVA, following the UK MILI "
    "design; a German trial now includes adolescents. A genotype-matched oligonucleotide "
    "depends on the child's own alleles, judged by rules we fixed in advance. Senolytics "
    "come last. Every row carries its main weakness, because a hypothesis without one isn't "
    "ready for an expert.",

    "A child with MVA may need cancer treatment at short notice, so every idea must also "
    "pass a six-point compatibility screen: no further weakening of the checkpoint, no "
    "CYP3A4 clash with vincristine, and no harm to kidney or bone marrow. We built an open "
    "pipeline that runs these checks over public data. On its first run it excluded every "
    "checkpoint-kinase inhibitor we fed it as a control, with no hand curation.",

    "Nothing in the rule depends on the gene. The same pipeline runs for TRIP13 and CEP57 "
    "from one line of configuration. The path is patient-derived cells first, then the "
    "five-assay rule, and only then a formal N-of-1 design.",

    "Veto first, then hypotheses for the experts. Thank you.",
]

SECTION_TITLES = [
    'The veto rule',
    'A dosage disease',
    'Nothing exists, and the databases point the wrong way',
    'Five hypotheses, each with its weakness',
    'Two filters, and a pipeline that applies them',
    'How it scales and where it goes next',
    'Close: veto first',
]

# Organisers' required acknowledgement (closing card only, not spoken).
ACK = ('This work was made possible through the Hackathon, organized by Sage Bionetworks in '
       'partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, '
       'Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and '
       'Anthropic. We are deeply grateful to the child and their family who generously '
       'contributed their data and their story to advance research into this rare disease. We '
       'acknowledge their trust in making this Hackathon possible.')

# Pronunciation hints. The words spoken are unchanged; these only tell the TTS how to say
# acronyms, gene symbols and a few technical words (misaki phonemes, American English).
KOKORO_LEXICON = {
    'MVA': 'ˌɛmvˌiˈA',
    'BUB1B': 'bˌʌb wˌʌn bˈi',
    'BUBR1': 'bˌʌb ˌɑɹ wˈʌn',
    'HDAC2/3': 'ˌAʧdˈæk tˈu ænd θɹˈi',
    'NAD+': 'ˌɛnˌAdˈi plˈʌs',
    'CYP3A4': 'sˈɪp θɹˈi ˈA fˈɔɹ',
    'TRIP13': 'tɹˈɪp θˌɜɹtˈin',
    'CEP57': 'sˈɛp fˈɪfti sˈɛvən',
    'ChEMBL': 'kˈɛmbəl',
    'PubMed': 'pˈʌbmˌɛd',
    'MILI': 'mˈɪli',
    'N-of-1': 'ˈɛn əv wˈʌn',
    'pseudokinase': 'sˌudOkˈInAs',
    'mitosis': 'mItˈOsɪs',
    'mitoses': 'mItˈOsiz',
    'acetylation': 'əsˌɛtᵊlˈAʃən',
    'acetylate': 'əsˈɛtᵊlˌAt',
    'Senolytics': 'sˌɛnəlˈɪtɪks',
    'Metformin': 'mɛtfˈɔɹmɪn',
}

# Plain-text respellings for the Piper fallback (Piper has no phoneme markup).
PIPER_RESPELL = {
    'MVA': 'M.V.A.', 'BUB1B': 'bub one B', 'BUBR1': 'bub R one',
    'HDAC2/3': 'H-dack two and three', 'NAD+': 'N.A.D. plus', 'CYP3A4': 'sip three A four',
    'TRIP13': 'trip thirteen', 'CEP57': 'sep fifty-seven', 'ChEMBL': 'kembl',
    'PubMed': 'Pub Med', 'MILI': 'Milly', 'N-of-1': 'N of one',
}

# ----------------------------------------------------------------------------- utilities
LOG = []


def log(msg):
    line = f'[{time.time() - T0:7.1f}s] {msg}'
    print(line, flush=True)
    LOG.append(line)


def run(cmd, check=True, timeout=None):
    shell = isinstance(cmd, str)
    r = subprocess.run(cmd, shell=shell, capture_output=True, text=True, timeout=timeout)
    if check and r.returncode != 0:
        shown = cmd if shell else ' '.join(map(str, cmd))
        raise RuntimeError(f'command failed ({r.returncode}): {shown[:600]}\n'
                           f'{(r.stdout or "")[-1500:]}\n{(r.stderr or "")[-3000:]}')
    return r


def pkg_version(name):
    try:
        import importlib.metadata as md
        return md.version(name)
    except Exception:
        return None


def mmss(t):
    return f'{int(t // 60)}:{t % 60:04.1f}'


def split_sentences(text):
    return [s.strip() for s in re.split(r'(?<=[.!?])\s+', text.strip()) if s.strip()]


def count_words(text):
    return len(re.findall(r"[A-Za-z0-9][A-Za-z0-9'+/\-]*", text))


def _lexicon_regex(keys):
    alt = '|'.join(re.escape(k) for k in sorted(keys, key=len, reverse=True))
    return re.compile(r'(?<![A-Za-z0-9])(' + alt + r')(?![A-Za-z0-9])')


KOKORO_RE = _lexicon_regex(KOKORO_LEXICON)
PIPER_RE = _lexicon_regex(PIPER_RESPELL)


def to_kokoro_input(sentence):
    return KOKORO_RE.sub(lambda m: f'[{m.group(1)}](/{KOKORO_LEXICON[m.group(1)]}/)', sentence)


def to_piper_input(sentence):
    return PIPER_RE.sub(lambda m: PIPER_RESPELL[m.group(1)], sentence)


# Heavy libraries are imported only after the TTS packages are installed, so that pip
# cannot swap a library that this process has already loaded.
np = Image = ImageDraw = ImageFont = RESAMPLE = None


def load_libs():
    global np, Image, ImageDraw, ImageFont, RESAMPLE
    import numpy as _np
    from PIL import Image as _Image, ImageDraw as _Draw, ImageFont as _Font
    np, Image, ImageDraw, ImageFont = _np, _Image, _Draw, _Font
    RESAMPLE = getattr(Image, 'Resampling', Image).LANCZOS


# ----------------------------------------------------------------------------- fonts
FONT_PREFS = ['Inter', 'Source Sans 3', 'Source Sans Pro', 'Open Sans', 'Noto Sans', 'Lato',
              'Roboto', 'Liberation Sans', 'DejaVu Sans']
SPECIAL_CHARS = '≈µ—“”’·§ '
NBSP = ' '  # keeps word pairs such as "EU CTIS" on one line
SYM_CHARS = '→'  # drawn with DejaVu Sans, whose arrow is larger and clearer
STYLE_KEYS = {'regular': 'regular', 'book': 'regular', 'normal': 'regular', 'roman': 'regular',
              'bold': 'bold', 'italic': 'italic', 'oblique': 'italic',
              'bolditalic': 'bolditalic', 'boldoblique': 'bolditalic'}


class Fonts:
    def __init__(self, name, paths, sym_paths):
        self.name, self.paths, self.sym_paths, self.cache = name, paths, sym_paths, {}

    def get(self, size, bold=False, italic=False, sym=False):
        if sym:
            path = self.sym_paths['bold' if bold else 'regular']
        else:
            style = ('bolditalic' if bold else 'italic') if italic else ('bold' if bold else 'regular')
            path = self.paths.get(style) or self.paths['bold' if bold else 'regular']
        px = max(1, int(round(size * SS)))
        key = (path, px)
        if key not in self.cache:
            self.cache[key] = ImageFont.truetype(path, px)
        return self.cache[key]


def pick_fonts():
    from fontTools.ttLib import TTFont
    required = set(string.ascii_letters + string.digits + " .,:;!?()/+-'&%") | set(SPECIAL_CHARS)
    fams = {}
    try:
        r = run(['fc-list', '--format', '%{family[0]}|%{style[0]}|%{file}\n'], check=False)
        for line in (r.stdout or '').splitlines():
            parts = line.split('|')
            if len(parts) != 3:
                continue
            fam, style, path = parts
            key = STYLE_KEYS.get(style.lower().replace(' ', ''))
            if key and path.lower().endswith(('.ttf', '.otf')):
                fams.setdefault(fam, {}).setdefault(key, path)
    except Exception as e:  # fc-list missing: fall back to matplotlib's DejaVu Sans
        log(f'fc-list unavailable: {e}')
    cands = [(n, fams[n]) for n in FONT_PREFS if n in fams]
    import matplotlib
    d = Path(matplotlib.get_data_path()) / 'fonts' / 'ttf'
    dejavu = {'regular': str(d / 'DejaVuSans.ttf'), 'bold': str(d / 'DejaVuSans-Bold.ttf'),
              'italic': str(d / 'DejaVuSans-Oblique.ttf'),
              'bolditalic': str(d / 'DejaVuSans-BoldOblique.ttf')}
    cands.append(('DejaVu Sans (matplotlib)', dict(dejavu)))
    for name, paths in cands:
        if 'regular' not in paths or 'bold' not in paths or paths['regular'] == paths['bold']:
            continue
        ok = True
        for k in ('regular', 'bold'):
            try:
                cmap = TTFont(paths[k], lazy=True).getBestCmap() or {}
            except Exception as e:
                log(f'font {name} {k}: unreadable ({e})')
                ok = False
                break
            missing = ''.join(sorted(c for c in required if ord(c) not in cmap))
            if missing:
                log(f'font {name} {k}: lacks {missing!r}, skipped')
                ok = False
                break
        if ok:
            paths = {k: v for k, v in paths.items() if Path(v).is_file()}
            log(f'font: {name} ({paths}); arrows from DejaVu Sans')
            return Fonts(name, paths, dejavu)
    raise RuntimeError('no usable font found')


# ----------------------------------------------------------------------------- drawing
C = {
    'bg': (11, 26, 44), 'panel': (19, 40, 64), 'panel2': (27, 53, 84), 'rule': (52, 80, 112),
    'text': (255, 255, 255), 'muted': (198, 210, 226), 'accent': (255, 193, 69),
    'teal': (88, 214, 200), 'warn': (255, 125, 108), 'ink': (11, 26, 44),
}
MARK_RE = re.compile(r'(\*\*.+?\*\*|\^\^.+?\^\^|\*[^*\s][^*]*?\*)')
TOKEN_RE = re.compile(r'[^\S ]+|(?:[^\s]| )+')  # a no-break space never splits


def parse(s, color='text', bold=False):
    """Tiny markup: **bold**, *italic*, ^^bold accent^^. Straight apostrophes become curly."""
    s = s.replace("'", '’')
    runs = []
    for part in MARK_RE.split(s):
        if not part:
            continue
        if part.startswith('**') and part.endswith('**') and len(part) > 4:
            run_ = (part[2:-2], True, False, color)
        elif part.startswith('^^') and part.endswith('^^') and len(part) > 4:
            run_ = (part[2:-2], True, False, 'accent')
        elif part.startswith('*') and part.endswith('*') and len(part) > 2:
            run_ = (part[1:-1], bold, True, color)
        else:
            run_ = (part, bold, False, color)
        for piece in re.split('([' + SYM_CHARS + '])', run_[0]):
            if piece:
                runs.append((piece,) + run_[1:] + (piece in SYM_CHARS,))
    return runs


def layout(runs, size, max_w, fonts):
    """Greedy word wrap. Words glued without a space (e.g. 'BUBR1**,') never split."""
    groups, pending_space = [], False
    for text, bold, italic, color, sym in runs:
        font = fonts.get(size, bold, italic, sym)
        for m in TOKEN_RE.finditer(text):
            t = m.group(0)
            if t.isspace():
                pending_space = True
                continue
            item = (t, font, C[color], font.getlength(t))
            if pending_space or not groups:
                groups.append({'space': pending_space and bool(groups), 'words': [item]})
            else:
                groups[-1]['words'].append(item)
            pending_space = False
    lines, cur, x = [], [], 0.0
    limit = max_w * SS
    for g in groups:
        gw = sum(it[3] for it in g['words'])
        sp = g['words'][0][1].getlength(' ') if (g['space'] and cur) else 0.0
        if cur and x + sp + gw > limit:
            lines.append(cur)
            cur, x, sp = [], 0.0, 0.0
        x += sp
        for t, font, col, w in g['words']:
            cur.append((x, t, font, col))
            x += w
    if cur:
        lines.append(cur)
    return lines


def line_width(line):
    if not line:
        return 0.0
    x, t, font, _ = line[-1]
    return x + font.getlength(t)


class Canvas:
    """All coordinates are in 1920x1080 units; drawing happens at SS x that size."""

    def __init__(self, fonts):
        self.f = fonts
        self.img = Image.new('RGB', (W * SS, H * SS), C['bg'])
        self.d = ImageDraw.Draw(self.img)
        self.bottom = 0.0

    def _track(self, y):
        self.bottom = max(self.bottom, y)

    def rect(self, x0, y0, x1, y1, fill=None, outline=None, width=0, r=0, track=True):
        box = [round(x0 * SS), round(y0 * SS), round(x1 * SS), round(y1 * SS)]
        kw = dict(fill=C[fill] if fill else None, outline=C[outline] if outline else None,
                  width=max(1, round(width * SS)) if outline else 0)
        if r > 0:
            self.d.rounded_rectangle(box, radius=round(r * SS), **kw)
        else:
            self.d.rectangle(box, **kw)
        if track:
            self._track(y1)

    def lines_for(self, s, size, max_w, color='text', bold=False):
        """Wrap, then narrow the measure as far as the line count allows (balanced lines,
        so no single word is left alone on the last line)."""
        def wrap(width):
            out = []
            for para in s.split('\n'):
                out.extend(layout(parse(para, color, bold), size, width, self.f) or [[]])
            return out

        lines = wrap(max_w)
        if len(lines) > 1:
            n, lo, hi = len(lines), max_w * 0.5, max_w
            for _ in range(10):
                mid = (lo + hi) / 2
                if len(wrap(mid)) <= n:
                    hi = mid
                else:
                    lo = mid
            lines = wrap(hi)
        return lines

    def block_height(self, n, size, lh, bold=False):
        asc, desc = self.f.get(size, bold).getmetrics()
        return (n - 1) * size * lh + (asc + desc) / SS

    def measure(self, s, size, max_w, bold=False, lh=1.25):
        lines = self.lines_for(s, size, max_w, 'text', bold)
        return self.block_height(len(lines), size, lh, bold), lines

    def text(self, s, x, y, size, max_w, color='text', bold=False, lh=1.25, align='left',
             track=True):
        lines = self.lines_for(s, size, max_w, color, bold)
        asc, _ = self.f.get(size, bold).getmetrics()
        for i, line in enumerate(lines):
            if not line:
                continue
            lw = line_width(line)
            if align == 'center':
                ox = x * SS + (max_w * SS - lw) / 2
            elif align == 'right':
                ox = x * SS + max_w * SS - lw
            else:
                ox = x * SS
            base = y * SS + asc + i * size * lh * SS
            for dx, t, font, col in line:
                self.d.text((ox + dx, base), t, font=font, fill=col, anchor='ls')
        h = self.block_height(len(lines), size, lh, bold)
        if track:
            self._track(y + h)
        return h

    def spaced(self, s, x, y, size, color, bold=True, tracking=0.12, right=False, track=True):
        font = self.f.get(size, bold)
        asc, desc = font.getmetrics()
        gap = tracking * size * SS
        widths = [font.getlength(ch) for ch in s]
        total = sum(widths) + gap * (len(s) - 1)
        ox = x * SS - total if right else x * SS
        for ch, w in zip(s, widths):
            self.d.text((ox, y * SS + asc), ch, font=font, fill=C[color], anchor='ls')
            ox += w + gap
        h = (asc + desc) / SS
        if track:
            self._track(y + h)
        return h

    def center_text(self, s, cx, cy, size, color, bold=True):
        self.d.text((cx * SS, cy * SS), s, font=self.f.get(size, bold), fill=C[color], anchor='mm')

    def box(self, x0, y0, x1, y1, s, size, fill='panel2', outline='accent', color='text',
            bold=True, r=18, width=3, pad=18, lh=1.16):
        self.rect(x0, y0, x1, y1, fill=fill, outline=outline, width=width, r=r)
        h, _ = self.measure(s, size, x1 - x0 - 2 * pad, bold=bold, lh=lh)
        ty = y0 + (y1 - y0 - h) / 2 + size * 0.04
        self.text(s, x0 + pad, ty, size, x1 - x0 - 2 * pad, color=color, bold=bold, lh=lh,
                  align='center')

    def arrow(self, x0, x1, y, color='accent', width=6, head=22):
        self.d.line([(x0 * SS, y * SS), ((x1 - head * 0.8) * SS, y * SS)], fill=C[color],
                    width=round(width * SS))
        self.d.polygon([(x1 * SS, y * SS), ((x1 - head) * SS, (y - head * 0.62) * SS),
                        ((x1 - head) * SS, (y + head * 0.62) * SS)], fill=C[color])

    def warn_icon(self, x, y, s):
        self.d.polygon([((x + s / 2) * SS, y * SS), ((x + s) * SS, (y + s * 0.9) * SS),
                        (x * SS, (y + s * 0.9) * SS)], fill=C['accent'])
        self.d.text(((x + s / 2) * SS, (y + s * 0.85) * SS), '!', font=self.f.get(s * 0.66, True),
                    fill=C['ink'], anchor='ms')

    def pill(self, s, x, y, size, fill='warn', color='ink', padx=30, pady=15, max_w=1700):
        inner = max_w - 2 * padx
        h, lines = self.measure(s, size, inner, bold=True, lh=1.2)
        w = max(line_width(l) for l in lines) / SS + 2 * padx
        hh = h + 2 * pady
        self.rect(x, y, x + w, y + hh, fill=fill, r=hh / 2 if len(lines) == 1 else 20)
        self.text(s, x + padx, y + pady, size, inner, color=color, bold=True, lh=1.2)
        return y + hh

    def bullets(self, items, x, y, w, size, gap=None, color='text', lh=1.25):
        gap = size * 0.5 if gap is None else gap
        indent = size * 0.9
        asc, _ = self.f.get(size).getmetrics()
        r = size * 0.16
        for it in items:
            cy = y + asc / SS - size * 0.33
            self.d.ellipse([round(x * SS), round((cy - r) * SS), round((x + 2 * r) * SS),
                            round((cy + r) * SS)], fill=C['accent'])
            h = self.text(it, x + indent, y, size, w - indent, color=color, lh=lh)
            y += h + gap
        return y - gap

    def title(self, s, y=100):
        size = None
        for cand in (68, 64, 60):
            _, lines = self.measure(s, cand, 1700, bold=True, lh=1.12)
            if len(lines) == 1:
                size = cand
                break
        size = size or 64
        h = self.text(s, 110, y, size, 1700, bold=True, lh=1.12)
        by = y + h + 20
        self.rect(110, by, 250, by + 8, fill='accent')
        return by + 8 + 44

    def chrome(self, counter=None, kicker=True):
        self.rect(0, 0, W, 6, fill='accent', track=False)
        if kicker:
            self.spaced(KICKER, 110, 44, 24, 'accent', track=False)
        if counter:
            self.spaced(counter, W - 110, 44, 24, 'muted', right=True, track=False)
        self.rect(110, 986, W - 110, 988, fill='rule', track=False)
        self.text(FOOTER, 110, 1004, 28, 1500, color='muted', track=False)


# ----------------------------------------------------------------------------- slides
def slide_title(cv, fs):
    cv.spaced('RARE DISEASE, REAL KID: THE MVA HACKATHON 2026', 110, 196, 30, 'accent',
              tracking=0.10)
    cv.spaced('TRACK 2 · DRUG REPURPOSING', 110, 248, 26, 'muted', tracking=0.12)
    h = cv.text('Veto first', 110, 318, 150 * fs, 1700, bold=True, lh=1.0)
    y = 318 + h + 30
    cv.rect(110, y, 290, y + 10, fill='accent')
    y += 10 + 44
    h = cv.text('Repurposing hypotheses for *BUB1B*-related mosaic variegated aneuploidy (MVA1)',
                110, y, 50 * fs, 1500, lh=1.22)
    y += h + 38
    cv.text('One veto rule · five hypotheses · two safety filters · one open pipeline',
            110, y, 34 * fs, 1700, color='muted')


def slide_1(cv, fs):
    y = cv.title(SECTION_TITLES[0])
    x0, x1, bar, pad = 110, W - 110, 14, 40
    quote = '“Any drug that shortens mitosis is out — whatever else it does.”'
    src = 'The five-assay decision rule (Track 2 report, §7)'
    inner = x1 - x0 - bar - 2 * pad
    qh, _ = cv.measure(quote, 56 * fs, inner, bold=True, lh=1.18)
    sh, _ = cv.measure(src, 30 * fs, inner)
    card_h = pad + qh + 20 + sh + pad
    cv.rect(x0, y, x1, y + card_h, fill='panel')
    cv.rect(x0, y, x0 + bar, y + card_h, fill='warn')
    cv.text(quote, x0 + bar + pad, y + pad, 56 * fs, inner, bold=True, lh=1.18)
    cv.text(src, x0 + bar + pad, y + pad + qh + 20, 30 * fs, inner, color='muted')
    y += card_h + 46 * fs
    cv.bullets([
        'Every drug in this pitch has to survive this rule **first**',
        "Shortens mitosis in the patient's own cells → out, however good its other results",
        'Comes straight from the biology of the disease',
        'Rules out several obvious drug classes before we start',
    ], 110, y, 1700, 44 * fs, gap=20 * fs)


def slide_2(cv, fs):
    y = cv.title(SECTION_TITLES[1])
    y = cv.bullets([
        'One *BUB1B* copy makes no protein; the other makes too little',
        'Adding BUBR1 back rescues patient cells',
        '**Goal: more BUBR1**, not a stronger kinase (BUBR1 is a pseudokinase)',
    ], 110, y, 1700, 44 * fs, gap=18 * fs) + 36 * fs
    h = cv.text('An acetylation cycle sets the BUBR1 level', 110, y, 32 * fs, 1700,
                color='accent', bold=True)
    y += h + 20 * fs
    n, gap = 3, 84
    bw = (1700 - (n - 1) * gap) / n
    bh = 148 * fs
    labels = ['PCAF acetylates K250:\nBUBR1 protected', 'HDAC2/3 remove\nthe mark',
              'APC/C destroys\nBUBR1']
    for i, lab in enumerate(labels):
        bx = 110 + i * (bw + gap)
        cv.box(bx, y, bx + bw, y + bh, lab, 34 * fs)
        if i < n - 1:
            cv.arrow(bx + bw + 14, bx + bw + gap - 14, y + bh / 2)
    y += bh + 40 * fs
    cv.pill('Acetylation-deficient BubR1 → shorter mitosis → MVA-like mice (Park 2013)',
            110, y, 34 * fs)


def slide_3(cv, fs):
    y = cv.title(SECTION_TITLES[2])
    tiles = [
        ('0', 'trials in 5 registries',
         f'ClinicalTrials.gov{NBSP}· ISRCTN{NBSP}· EU{NBSP}CTIS{NBSP}· ANZCTR{NBSP}· WHO{NBSP}ICTRP'),
        ('0', 'disease-modifying reports', 'PubMed'),
        ('≈10 µM', 'best BUB1B ligand', 'ChEMBL'),
    ]
    gap, pad = 34, 32
    tw = (1700 - 2 * gap) / 3
    inner = tw - 2 * pad
    big, lab, small = 100 * fs, 34 * fs, 26 * fs
    hb = max(cv.measure(b, big, inner, bold=True, lh=1.0)[0] for b, _, _ in tiles)
    hl = max(cv.measure(l, lab, inner, bold=True)[0] for _, l, _ in tiles)
    hs = max(cv.measure(s, small, inner)[0] for _, _, s in tiles)
    th = pad + hb + 18 + hl + 12 + hs + pad
    for i, (b, l, s) in enumerate(tiles):
        tx = 110 + i * (tw + gap)
        cv.rect(tx, y, tx + tw, y + th, fill='panel', r=18)
        cv.text(b, tx + pad, y + pad, big, inner, color='accent', bold=True, lh=1.0)
        cv.text(l, tx + pad, y + pad + hb + 18, lab, inner, bold=True)
        cv.text(s, tx + pad, y + pad + hb + 18 + hl + 12, small, inner, color='muted')
    y += th + 44 * fs
    h = cv.text('Open Targets: the only positive tractability signals are traps', 110, y,
                32 * fs, 1700, color='accent', bold=True)
    y += h + 24 * fs
    isz = 44 * fs
    for t in ['“Druggable family” is a kinase heuristic, but BUBR1 is a pseudokinase',
              '“Degrader-tractable” is the harmful direction for a dosage disease']:
        cv.warn_icon(110, y + 2 * fs, isz)
        h = cv.text(t, 110 + isz + 26, y, 40 * fs, 1700 - isz - 26)
        y += h + 22 * fs


def slide_4(cv, fs):
    y = cv.title(SECTION_TITLES[3])
    h = cv.text('Ranked by how cheaply the first experiment could disprove them', 110, y - 10,
                36 * fs, 1700, color='muted')
    y += h + 22 * fs - 10
    rows = [
        ('H1', 'NAD+ precursors via SIRT2', 'SIRT2 also removes the protective mark'),
        ('H2', 'HDAC2/3 inhibition', 'The class can cause aneuploidy and low platelets'),
        ('H3', 'Metformin inside surveillance (UK MILI design)',
         "May select against the child's own mosaic tissue"),
        ('H4', 'Genotype-matched oligonucleotide', "Depends on the child's alleles"),
        ('H5', 'Senolytics', 'Senescence may be a brake on tumours'),
    ]
    xb, xh, xw, xr = 128, 230, 1040, W - 128
    wh, ww = xw - 40 - xh, xr - xw
    cv.spaced('HYPOTHESIS', xh, y, 24 * fs, 'accent')
    cv.spaced('MAIN WEAKNESS', xw, y, 24 * fs, 'accent')
    y += 24 * fs * 1.25 + 14 * fs
    size, vpad, rgap, bs = 36 * fs, 14 * fs, 10 * fs, 58 * fs
    for i, (tag, hyp, weak) in enumerate(rows):
        hh, _ = cv.measure(hyp, size, wh, bold=True)
        hw, _ = cv.measure(weak, size, ww)
        rh = max(hh, hw, bs) + 2 * vpad
        cv.rect(110, y, W - 110, y + rh, fill='panel' if i % 2 == 0 else 'panel2', r=14)
        by = y + (rh - bs) / 2
        cv.rect(xb, by, xb + bs, by + bs, fill='accent', r=12)
        cv.center_text(tag, xb + bs / 2, by + bs / 2, 28 * fs, 'ink')
        cv.text(hyp, xh, y + (rh - hh) / 2, size, wh, bold=True)
        cv.text(weak, xw, y + (rh - hw) / 2, size, ww)
        y += rh + rgap


def slide_5(cv, fs):
    y = cv.title(SECTION_TITLES[4])
    y = cv.bullets(['A child with MVA may need cancer treatment at short notice'], 110, y, 1700,
                   44 * fs) + 36 * fs
    pad, gapx = 34, 40
    lx0, lx1 = 110, 910
    rx0, rx1 = lx1 + gapx, W - 110
    lab1 = 'Six-point chemotherapy-\ncompatibility screen'
    chips = ['Checkpoint', 'CYP3A4', 'Kidney', 'Marrow', 'Genotoxicity', 'Quick to stop']
    linner = lx1 - lx0 - 2 * pad
    hl1, _ = cv.measure(lab1, 32 * fs, linner, bold=True)
    cg = 22
    cw = (linner - cg) / 2
    ch = 76 * fs
    left_h = pad + hl1 + 22 * fs + 3 * ch + 2 * cg + pad
    lab2 = 'Open pipeline, first run on public data'
    rinner = rx1 - rx0 - 2 * pad
    hl2, _ = cv.measure(lab2, 32 * fs, rinner, bold=True)
    bars = [('28 of 28', 'checkpoint-kinase controls excluded automatically', 1.0),
            ('11 of 12', 'HDAC2/3 candidates excluded, with no hand curation', 11 / 12)]
    bar_h = 20 * fs
    rows_h = [(cv.measure(num, 48 * fs, rinner, bold=True)[0], cv.measure(desc, 30 * fs, rinner)[0])
              for num, desc, _ in bars]
    right_h = (pad + hl2 + 22 * fs + sum(hn + 8 + hd + 14 + bar_h for hn, hd in rows_h)
               + 26 * fs * (len(bars) - 1) + pad)
    ph = max(left_h, right_h)
    cv.rect(lx0, y, lx1, y + ph, fill='panel', r=18)
    cv.rect(rx0, y, rx1, y + ph, fill='panel', r=18)
    yy = y + pad
    cv.text(lab1, lx0 + pad, yy, 32 * fs, linner, color='accent', bold=True)
    yy += hl1 + 22 * fs
    for i, chip in enumerate(chips):
        cx = lx0 + pad + (i % 2) * (cw + cg)
        cy = yy + (i // 2) * (ch + cg)
        cv.box(cx, cy, cx + cw, cy + ch, chip, 32 * fs, fill='panel2', outline='teal', r=ch / 2)
    yy = y + pad
    cv.text(lab2, rx0 + pad, yy, 32 * fs, rinner, color='accent', bold=True)
    yy += hl2 + 22 * fs
    for (num, desc, frac), (hn, hd) in zip(bars, rows_h):
        cv.text(num, rx0 + pad, yy, 48 * fs, rinner, bold=True)
        yy += hn + 8
        cv.text(desc, rx0 + pad, yy, 30 * fs, rinner, color='muted')
        yy += hd + 14
        cv.rect(rx0 + pad, yy, rx0 + pad + rinner, yy + bar_h, fill='rule', r=bar_h / 2)
        cv.rect(rx0 + pad, yy, rx0 + pad + rinner * frac, yy + bar_h, fill='teal', r=bar_h / 2)
        yy += bar_h + 26 * fs


def slide_6(cv, fs):
    y = cv.title(SECTION_TITLES[5])
    y = cv.bullets(['Nothing in the rule depends on the gene'], 110, y, 1700, 44 * fs) + 40 * fs
    h = cv.text('Three genes, one config file', 110, y, 32 * fs, 1700, color='accent', bold=True)
    y += h + 20 * fs
    cw, ch, cg = 300, 92 * fs, 30
    for i, g in enumerate(['*BUB1B*', '*TRIP13*', '*CEP57*']):
        cx = 110 + i * (cw + cg)
        cv.box(cx, y, cx + cw, y + ch, g, 42 * fs, fill='panel2', outline='teal', r=ch / 2)
    y += ch + 46 * fs
    h = cv.text('Where it goes next', 110, y, 32 * fs, 1700, color='accent', bold=True)
    y += h + 20 * fs
    steps = ['Patient-derived\ncells', 'Five-assay\nrule', 'N-of-1 design\n(CENT)', 'Registry']
    n, gap = len(steps), 76
    bw = (1700 - (n - 1) * gap) / n
    bh = 132 * fs
    for i, s in enumerate(steps):
        bx = 110 + i * (bw + gap)
        cv.box(bx, y, bx + bw, y + bh, s, 34 * fs)
        if i < n - 1:
            cv.arrow(bx + bw + 14, bx + bw + gap - 14, y + bh / 2)
    y += bh + 44 * fs
    cv.bullets(['Adding a gene takes one line of configuration'], 110, y, 1700, 40 * fs)


def slide_7(cv, fs):
    h = cv.text('Veto first.', 110, 112, 104 * fs, 1700, color='accent', bold=True, lh=1.0)
    y = 112 + h + 16
    h = cv.text('Hypotheses, not prescriptions.', 110, y, 60 * fs, 1700, bold=True)
    y += h + 28
    if REPO_URL:
        h = cv.text('Code and reports: ' + REPO_URL, 110, y, 32 * fs, 1700, color='teal')
        y += h + 16
    y += 26
    pad, bar = 38, 12
    inner = 1700 - bar - 2 * pad
    lab_h = 24 * fs * 1.2
    ah, _ = cv.measure(ACK, 34 * fs, inner, lh=1.4)
    card_h = pad + lab_h + 18 + ah + pad
    cv.rect(110, y, W - 110, y + card_h, fill='panel')
    cv.rect(110, y, 110 + bar, y + card_h, fill='accent')
    cv.spaced('ACKNOWLEDGEMENT', 110 + bar + pad, y + pad, 24 * fs, 'accent')
    cv.text(ACK, 110 + bar + pad, y + pad + lab_h + 18, 34 * fs, inner, lh=1.4)


SLIDES = [  # key, title, renderer, counter, show kicker, largest body font scale
    ('title', 'Veto first (title card)', slide_title, None, False, 1.0),
    ('s1', SECTION_TITLES[0], slide_1, '1 / 7', True, 1.16),
    ('s2', SECTION_TITLES[1], slide_2, '2 / 7', True, 1.16),
    ('s3', SECTION_TITLES[2], slide_3, '3 / 7', True, 1.16),
    ('s4', SECTION_TITLES[3], slide_4, '4 / 7', True, 1.16),
    ('s5', SECTION_TITLES[4], slide_5, '5 / 7', True, 1.16),
    ('s6', SECTION_TITLES[5], slide_6, '6 / 7', True, 1.16),
    ('s7', SECTION_TITLES[6], slide_7, '7 / 7', True, 1.16),
]


def render_slides(fonts):
    (OUT / 'slides').mkdir(parents=True, exist_ok=True)
    info = []
    for idx, (key, title, fn, counter, kicker, max_fs) in enumerate(SLIDES):
        for fs in [s for s in FONT_SCALES if s <= max_fs]:
            cv = Canvas(fonts)
            fn(cv, fs)
            if cv.bottom <= CONTENT_LIMIT:
                break
        overflow = cv.bottom > CONTENT_LIMIT
        cv.chrome(counter, kicker)
        path = OUT / 'slides' / f'slide_{idx:02d}_{key}.png'
        cv.img.resize((W, H), RESAMPLE).save(path)
        info.append(dict(index=idx, key=key, title=title, png=f'slides/{path.name}',
                         font_scale=fs, content_bottom_px=round(cv.bottom, 1), overflow=overflow))
        log(f'slide {idx} {key}: font scale {fs}, content bottom {cv.bottom:.0f}px'
            + (' OVERFLOW' if overflow else ''))
    return info


# ----------------------------------------------------------------------------- TTS
def pip_install(args, py=None, timeout=1800):
    return run([py or sys.executable, '-m', 'pip', 'install', '-q', '--disable-pip-version-check',
                *args], timeout=timeout)


def apt_install_espeak():
    r = run('(apt-get -qq update && DEBIAN_FRONTEND=noninteractive apt-get -qq install -y '
            'espeak-ng) > /tmp/mva_apt.log 2>&1', check=False, timeout=900)
    log(f'apt-get install espeak-ng: exit {r.returncode}')
    return r.returncode == 0


def setup_kokoro_host():
    """Path A: Kokoro in the kernel's own Python. Every kokoro release so far declares
    Requires-Python <3.13, so on 3.13+ the code is installed without that check and its
    runtime dependencies are added by hand (spacy-curated-transformers is skipped: only the
    transformer spaCy model needs it, and misaki uses the small model)."""
    if sys.version_info < (3, 13):
        pip_install(['kokoro>=0.9.4', 'misaki[en]>=0.9.4', 'soundfile'])
    else:
        pip_install(['--ignore-requires-python', '--no-deps', 'kokoro>=0.9.4', 'misaki>=0.9.4'],
                    timeout=900)
        pip_install(['loguru', 'num2words', 'phonemizer-fork', 'espeakng-loader', 'spacy',
                     'soundfile', 'addict', 'regex'], timeout=900)
    run([sys.executable, '-m', 'spacy', 'download', 'en_core_web_sm'], timeout=900)
    log('Kokoro packages installed in the kernel Python')
    return sys.executable


def setup_kokoro_venv():
    """Path B: a uv-managed Python 3.12 venv with the versions kokoro 0.9.4 was released
    against (CPU-only torch)."""
    venv = Path('/tmp/kokoro-venv')
    pip_install(['uv'], timeout=600)
    uv = [sys.executable, '-m', 'uv']
    run(uv + ['venv', '--seed', '--python', '3.12', str(venv)], timeout=900)
    py = str(venv / 'bin' / 'python')
    run(uv + ['pip', 'install', '--python', py, 'torch==2.6.0',
              '--index-url', 'https://download.pytorch.org/whl/cpu'], timeout=1800)
    run(uv + ['pip', 'install', '--python', py, 'kokoro==0.9.4', 'misaki[en]==0.9.4',
              'transformers<5', 'soundfile', 'numpy'], timeout=1800)
    run([py, '-m', 'spacy', 'download', 'en_core_web_sm'], timeout=900)
    log('Kokoro packages installed in a Python 3.12 venv')
    return py


# Runs in its own Python process, so Kokoro's dependencies never mix with this kernel's.
KOKORO_WORKER = r'''
import json, os, re, sys, time
t0 = time.time()
job = json.load(open(sys.argv[1], encoding="utf-8"))
import numpy as np
import torch
torch.set_num_threads(max(1, os.cpu_count() or 1))
from kokoro import KPipeline
pipe = KPipeline(lang_code=job["lang"], repo_id="hexgrad/Kokoro-82M")
os.makedirs(job["out_dir"], exist_ok=True)
results = []
for item in job["items"]:
    chunks, phon, unknown = [], [], []
    for res in pipe(item["input"], voice=job["voice"], speed=job["speed"], split_pattern=None):
        if isinstance(res, tuple):
            ps, audio, toks = res[1], res[2], []
        else:
            ps, audio, toks = res.phonemes, res.audio, (getattr(res, "tokens", None) or [])
        for tk in toks:
            ph = getattr(tk, "phonemes", None)
            head = getattr(getattr(tk, "_", None), "is_head", True)
            if not ph and head and re.search(r"[A-Za-z0-9]", getattr(tk, "text", "") or ""):
                unknown.append(tk.text)
        if audio is None:
            continue
        a = audio.detach().cpu().numpy() if hasattr(audio, "detach") else np.asarray(audio)
        chunks.append(a.astype(np.float32).reshape(-1))
        phon.append(ps)
    if not chunks:
        sys.exit("Kokoro returned no audio for item " + item["id"])
    path = os.path.join(job["out_dir"], item["id"] + ".npy")
    np.save(path, np.concatenate(chunks))
    results.append({"id": item["id"], "npy": path, "phonemes": " | ".join(phon),
                    "unknown": unknown})


def ver(p):
    try:
        import importlib.metadata as md
        return md.version(p)
    except Exception:
        return None


out = {"results": results, "python": sys.version.split()[0], "seconds": round(time.time() - t0, 1),
       "versions": {p: ver(p) for p in ("kokoro", "misaki", "torch", "transformers", "spacy")}}
with open(os.path.join(job["out_dir"], "result.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False)
'''


class KokoroEngine:
    """Kokoro-82M, run through KOKORO_WORKER with the given Python interpreter."""
    name = 'Kokoro-82M'

    def __init__(self, python, how):
        self.python, self.how = python, how
        self.sr, self.voice, self.versions, self.calls = 24000, VOICE, {}, 0
        self.worker = WORK / 'kokoro_worker.py'
        self.worker.write_text(KOKORO_WORKER, encoding='utf-8')

    def prepare(self, sentence):
        return to_kokoro_input(sentence)

    def synth_batch(self, sentences, speed):
        self.calls += 1
        out_dir = WORK / f'kokoro_{self.calls}'
        job = WORK / f'kokoro_job_{self.calls}.json'
        with open(job, 'w', encoding='utf-8') as f:
            json.dump({'lang': LANG, 'voice': self.voice, 'speed': speed, 'out_dir': str(out_dir),
                       'items': [{'id': f'{i:03d}', 'input': self.prepare(s)}
                                 for i, s in enumerate(sentences)]}, f, ensure_ascii=False)
        run([self.python, '-u', str(self.worker), str(job)], timeout=3600)
        with open(out_dir / 'result.json', encoding='utf-8') as f:
            res = json.load(f)
        self.versions = dict(res.get('versions') or {}, python=res.get('python'),
                             environment=self.how)
        log(f'Kokoro worker ({self.how}): {len(sentences)} sentences in {res.get("seconds")}s')
        return [(np.load(x['npy']).astype(np.float32).reshape(-1), x['phonemes'], x['unknown'])
                for x in res['results']]


class PiperEngine:
    name = 'Piper'

    def __init__(self):
        run([sys.executable, '-m', 'pip', 'install', '-q', '--disable-pip-version-check',
             'piper-tts'], timeout=1800)
        from huggingface_hub import hf_hub_download
        base = 'en/en_US/lessac/high/en_US-lessac-high'
        onnx = hf_hub_download('rhasspy/piper-voices', base + '.onnx')
        hf_hub_download('rhasspy/piper-voices', base + '.onnx.json')
        from piper import PiperVoice
        self.v = PiperVoice.load(onnx)
        self.sr = int(self.v.config.sample_rate)
        self.voice = PIPER_VOICE
        self.versions = {p: pkg_version(p) for p in ('piper-tts', 'onnxruntime')}
        log(f'Piper ready: {self.versions}')

    def prepare(self, sentence):
        return to_piper_input(sentence)

    def synth(self, sentence, speed):
        path = WORK / 'piper_tmp.wav'
        text = self.prepare(sentence)
        with wave.open(str(path), 'wb') as wf:
            if hasattr(self.v, 'synthesize_wav'):  # piper-tts >= 1.3
                try:
                    from piper import SynthesisConfig
                    self.v.synthesize_wav(text, wf, syn_config=SynthesisConfig(length_scale=1.0 / speed))
                except ImportError:
                    self.v.synthesize_wav(text, wf)
            else:                                  # piper-tts 1.2
                self.v.synthesize(text, wf, length_scale=1.0 / speed)
        with wave.open(str(path), 'rb') as wf:
            data = wf.readframes(wf.getnframes())
        return np.frombuffer(data, dtype=np.int16).astype(np.float32) / 32768.0, text, []

    def synth_batch(self, sentences, speed):
        return [self.synth(s, speed) for s in sentences]


def smoke_test(eng):
    (a, ph, unk), = eng.synth_batch(['MVA type 1 is a dosage disease.'], 1.0)
    if a.size < eng.sr * 0.5:
        raise RuntimeError(f'smoke test returned only {a.size} samples')
    log(f'{eng.name} smoke test: {a.size / eng.sr:.2f}s of audio, unknown words {unk}')


def make_engine(host_python, host_error):
    """Kokoro in the kernel Python (A), else Kokoro in a Python 3.12 venv (B), else Piper."""
    errors, attempts = [], []
    if host_python:
        attempts.append((f'kernel Python {platform.python_version()}', lambda: host_python))
    else:
        errors.append(f'kokoro (kernel Python): {host_error}')
    attempts.append(('uv venv, Python 3.12', setup_kokoro_venv))
    for how, setup in attempts:
        try:
            eng = KokoroEngine(setup(), how)
            smoke_test(eng)
            return eng, errors
        except Exception as e:
            errors.append(f'kokoro ({how}): {type(e).__name__}: {str(e)[-1500:]}')
            log(f'Kokoro via {how} failed:\n{traceback.format_exc()[-4000:]}')
    try:
        eng = PiperEngine()
        smoke_test(eng)
        return eng, errors
    except Exception as e:
        errors.append(f'piper: {type(e).__name__}: {str(e)[-1500:]}')
        log('Piper failed:\n' + traceback.format_exc()[-4000:])
    raise RuntimeError('no TTS engine worked: ' + ' | '.join(errors))


def trim_silence(a, sr, thresh_db=-45.0, pad_start=0.05, pad_end=0.09):
    win = max(1, int(0.01 * sr))
    n = len(a) // win
    if n < 3:
        return a
    rms = np.sqrt((a[:n * win].reshape(n, win) ** 2).mean(axis=1) + 1e-12)
    db = 20 * np.log10(rms / rms.max())
    idx = np.where(db > thresh_db)[0]
    if idx.size == 0:
        return a
    s = max(0, idx[0] * win - int(pad_start * sr))
    e = min(len(a), (idx[-1] + 1) * win + int(pad_end * sr))
    return a[s:e]


def synth_all(engine, speed):
    flat = [(si, s) for si, text in enumerate(NARRATION) for s in split_sentences(text)]
    outs = engine.synth_batch([s for _, s in flat], speed)
    secs = [[] for _ in NARRATION]
    for (si, s), (a, ph, unk) in zip(flat, outs):
        secs[si].append(dict(text=s, tts_input=engine.prepare(s), phonemes=ph, unknown=unk,
                             audio=trim_silence(a, engine.sr)))
    for si, sents in enumerate(secs):
        log(f'  section {si + 1}: {len(sents)} sentences, '
            f'{sum(len(x["audio"]) for x in sents) / engine.sr:.1f}s speech')
    return secs


def build_timeline(secs, sr, tm):
    def q(x):  # whole video frames
        return math.ceil(x * FPS - 1e-6) / FPS

    slides = [dict(index=0, start=0.0, dur=q(tm['title_hold']), audio=None, a_len=0.0, spans=[])]
    t = slides[0]['dur']
    gap = np.zeros(int(round(tm['sentence_pause'] * sr)), np.float32)
    for k, sents in enumerate(secs, start=1):
        parts, spans, off = [], [], 0
        for i, s in enumerate(sents):
            if i:
                parts.append(gap)
                off += len(gap)
            spans.append((off / sr, (off + len(s['audio'])) / sr))
            parts.append(s['audio'])
            off += len(s['audio'])
        audio = np.concatenate(parts)
        a_len = len(audio) / sr
        end_pad = tm['close_hold'] if k == len(secs) else tm['tail']
        d = q(tm['lead'] + a_len + end_pad)
        slides.append(dict(index=k, start=t, dur=d, audio=audio, a_len=a_len, spans=spans))
        t += d
    return slides, t


def plan_audio(engine):
    tm = dict(TIMING)
    speed = 1.0
    secs = synth_all(engine, speed)
    slides, total = build_timeline(secs, engine.sr, tm)
    log(f'pass 1 at speed {speed}: total {total:.2f}s')
    if total > MAX_TOTAL or total < MIN_TOTAL:
        speech = sum(len(s['audio']) for sec in secs for s in sec) / engine.sr
        other = total - speech
        goal = (MAX_TOTAL - 3.0) if total > MAX_TOTAL else TARGET_TOTAL
        new = min(SPEED_MAX, max(SPEED_MIN, speech / max(1.0, goal - other)))
        if abs(new - speed) >= 0.01:
            speed = round(new, 3)
            secs = synth_all(engine, speed)
            slides, total = build_timeline(secs, engine.sr, tm)
            log(f'pass 2 at speed {speed}: total {total:.2f}s')
    for key, val in [('close_hold', 10.0), ('tail', 0.55), ('title_hold', 3.0),
                     ('close_hold', 9.0), ('sentence_pause', 0.24), ('lead', 0.4)]:
        if total <= MAX_TOTAL:
            break
        tm[key] = val
        slides, total = build_timeline(secs, engine.sr, tm)
        log(f'shortened {key} to {val}: total {total:.2f}s')
    for key, val in [('close_hold', 14.0), ('tail', 0.85), ('sentence_pause', 0.36),
                     ('close_hold', 15.0), ('tail', 1.0)]:
        if total >= MIN_TOTAL:
            break
        tm[key] = val
        slides, total = build_timeline(secs, engine.sr, tm)
        log(f'lengthened {key} to {val}: total {total:.2f}s')
    if total >= HARD_MAX:
        raise RuntimeError(f'video would run {total:.1f}s even at speed {speed}; not cutting words')
    unknown = sorted({u for sec in secs for s in sec for u in s['unknown']})
    if unknown:
        log(f'WARNING: the TTS had no pronunciation for (and dropped): {unknown}')
    return secs, slides, total, speed, tm, unknown


def narration_track(slides, total, sr, lead):
    n = int(round(total * sr))
    track = np.zeros(n, np.float32)
    for sl in slides[1:]:
        off = int(round((sl['start'] + lead) * sr))
        a = sl['audio'][:max(0, n - off)]
        track[off:off + len(a)] = a
    peak = float(np.max(np.abs(track))) if n else 0.0
    if peak < 1e-4:
        raise RuntimeError('narration track is silent')
    return track * (0.89 / peak)  # peak about -1 dBFS


def write_wav(path, audio, sr):
    pcm = (np.clip(audio, -1.0, 1.0) * 32767.0).astype('<i2')
    with wave.open(str(path), 'wb') as wf:
        wf.setnchannels(1)
        wf.setsampwidth(2)
        wf.setframerate(sr)
        wf.writeframes(pcm.tobytes())


# ----------------------------------------------------------------------------- video
def find_ffmpeg():
    ff = shutil.which('ffmpeg')
    if not ff:
        run('DEBIAN_FRONTEND=noninteractive apt-get -qq install -y ffmpeg > /tmp/mva_apt_ff.log 2>&1',
            check=False, timeout=900)
        ff = shutil.which('ffmpeg')
    if not ff:
        try:
            import imageio_ffmpeg
            ff = imageio_ffmpeg.get_ffmpeg_exe()
        except Exception:
            ff = None
    if not ff:
        raise RuntimeError('ffmpeg not found')
    fp = shutil.which('ffprobe')
    ver = (run([ff, '-hide_banner', '-version'], check=False).stdout or '').splitlines()[:1]
    log(f'ffmpeg: {ff} {ver}; ffprobe: {fp}')
    if 'libx264' not in (run([ff, '-hide_banner', '-encoders'], check=False).stdout or ''):
        raise RuntimeError('this ffmpeg has no libx264 encoder')
    return ff, fp


def measure_loudness(ff, path):
    """EBU R128 integrated loudness, true peak and range (ffmpeg loudnorm, analysis only)."""
    txt = run([ff, '-hide_banner', '-nostats', '-i', str(path), '-af',
               'loudnorm=I=-16:TP=-1.5:LRA=11:print_format=json', '-f', 'null', '-'],
              check=False).stderr or ''
    m = re.search(r'\{[^{}]*"input_i"[^{}]*\}', txt, re.S)
    if not m:
        return {}
    d = json.loads(m.group(0))
    out = {}
    for k in ('input_i', 'input_tp', 'input_lra'):
        try:
            out[k] = float(d[k])
        except (KeyError, TypeError, ValueError):
            out[k] = None
    return out


def encode(ff, pngs, durs, wav, total, out, gain_db):
    cmd = [ff, '-y', '-hide_banner', '-loglevel', 'error']
    for png, d in zip(pngs, durs):
        cmd += ['-loop', '1', '-framerate', str(FPS), '-t', f'{d:.6f}', '-i', str(png)]
    cmd += ['-i', str(wav)]
    n = len(pngs)
    fo = max(0.0, total - 1.0)
    vf = (''.join(f'[{i}:v]' for i in range(n)) + f'concat=n={n}:v=1:a=0,'
          'scale=out_color_matrix=bt709:out_range=tv,format=yuv420p,'
          f'fade=t=in:st=0:d=0.6,fade=t=out:st={fo:.3f}:d=1.0[v]')
    # gain toward the loudness target, a look-ahead limiter at -1.5 dBFS, and the mono voice
    # copied at full level to both channels (ffmpeg's default upmix would drop it by 3 dB)
    af = (f'[{n}:a]aresample=48000,volume={gain_db:.2f}dB,'
          'alimiter=limit=0.84:attack=5:release=60:level=0:latency=1,'
          f'pan=stereo|c0=c0|c1=c0,afade=t=out:st={fo:.3f}:d=1.0[a]')
    cmd += ['-filter_complex', vf + ';' + af, '-map', '[v]', '-map', '[a]',
            '-c:v', 'libx264', '-preset', 'medium', '-tune', 'stillimage', '-crf', '18',
            '-profile:v', 'high', '-level:v', '4.1', '-pix_fmt', 'yuv420p', '-r', str(FPS),
            '-g', str(2 * FPS), '-colorspace', 'bt709', '-color_primaries', 'bt709',
            '-color_trc', 'bt709', '-color_range', 'tv',
            '-c:a', 'aac', '-b:a', '192k', '-ar', '48000', '-ac', '2',
            '-t', f'{total:.3f}', '-movflags', '+faststart', str(out)]
    run(cmd, timeout=3600)


def probe(ff, fp, path):
    if fp:
        r = run([fp, '-v', 'error', '-show_entries',
                 'format=duration,size,bit_rate:stream=index,codec_type,codec_name,profile,width,'
                 'height,pix_fmt,r_frame_rate,sample_rate,channels,duration',
                 '-of', 'json', str(path)])
        return json.loads(r.stdout)
    txt = run([ff, '-hide_banner', '-i', str(path)], check=False).stderr or ''
    m = re.search(r'Duration: (\d+):(\d+):([\d.]+)', txt)
    dur = int(m.group(1)) * 3600 + int(m.group(2)) * 60 + float(m.group(3)) if m else None
    return {'format': {'duration': dur, 'size': path.stat().st_size}, 'streams': [], 'raw': txt[-1500:]}


def volume_stats(ff, path):
    txt = run([ff, '-hide_banner', '-nostats', '-i', str(path), '-map', '0:a:0', '-af',
               'volumedetect', '-f', 'null', '-'], check=False).stderr or ''
    out = {}
    for key in ('mean_volume', 'max_volume'):
        m = re.search(key + r':\s*(-?[\d.]+) dB', txt)
        out[key + '_db'] = float(m.group(1)) if m else None
    return out


def grab_frame(ff, video, t, out_png, scale=None):
    cmd = [ff, '-y', '-v', 'error', '-ss', f'{t:.3f}', '-i', str(video), '-frames:v', '1']
    if scale:
        cmd += ['-vf', f'scale={scale}']
    run(cmd + [str(out_png)])


def contact_sheet(ff, video, slides, out_png):
    tw, th, gap, cols = 640, 360, 10, 4
    rows = math.ceil(len(slides) / cols)
    sheet = Image.new('RGB', (cols * tw + (cols + 1) * gap, rows * th + (rows + 1) * gap), (0, 0, 0))
    for i, sl in enumerate(slides):
        tmp = WORK / f'cs_{i}.png'
        grab_frame(ff, video, sl['start'] + sl['dur'] / 2, tmp, scale=f'{tw}:{th}')
        sheet.paste(Image.open(tmp).convert('RGB'),
                    (gap + (i % cols) * (tw + gap), gap + (i // cols) * (th + gap)))
    sheet.save(out_png)


def srt_ts(t):
    ms = int(round(t * 1000))
    h, ms = divmod(ms, 3600000)
    m, ms = divmod(ms, 60000)
    s, ms = divmod(ms, 1000)
    return f'{h:02d}:{m:02d}:{s:02d},{ms:03d}'


def caption_chunks(sentence, max_len=84):
    if len(sentence) <= max_len:
        return [sentence]
    chunks, cur = [], ''
    for p in re.split(r'(?<=[,;:])\s+', sentence):
        if cur and len(cur) + 1 + len(p) > max_len:
            chunks.append(cur)
            cur = p
        else:
            cur = (cur + ' ' + p).strip()
    if cur:
        chunks.append(cur)
    out = []
    for c in chunks:
        while len(c) > max_len:
            cut = c.rfind(' ', 0, max_len)
            if cut <= 0:
                break
            out.append(c[:cut])
            c = c[cut + 1:]
        out.append(c)
    return out


def write_captions(slides, secs, lead, path):
    cues = []
    for sl, sents in zip(slides[1:], secs):
        base = sl['start'] + lead
        for (s0, s1), s in zip(sl['spans'], sents):
            chunks = caption_chunks(s['text'])
            n_chars = sum(len(c) for c in chunks)
            acc = 0
            for c in chunks:
                c0 = base + s0 + (s1 - s0) * acc / n_chars
                acc += len(c)
                c1 = base + s0 + (s1 - s0) * acc / n_chars
                cues.append((c0, c1, c))
    with open(path, 'w', encoding='utf-8') as f:
        for i, (a, b, c) in enumerate(cues, 1):
            f.write(f'{i}\n{srt_ts(a)} --> {srt_ts(b)}\n' + '\n'.join(textwrap.wrap(c, 46)) + '\n\n')
    return len(cues)


# ----------------------------------------------------------------------------- main
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    summary = {'status': 'running', 'kernel': 'santanubanerjee9/mva-pitch-video',
               'generated_utc': datetime.now(timezone.utc).isoformat(timespec='seconds')}
    try:
        log(f'python {sys.version.split()[0]}, {os.cpu_count()} CPUs')
        apt_install_espeak()
        host_python, host_error = None, None
        try:                                          # before numpy/PIL are imported here
            host_python = setup_kokoro_host()
        except Exception as e:
            host_error = f'{type(e).__name__}: {str(e)[-1500:]}'
            log('Kokoro install in the kernel Python failed: ' + host_error)
        load_libs()
        ff, fp = find_ffmpeg()
        fonts = pick_fonts()
        summary['font'] = fonts.name
        slides_info = render_slides(fonts)

        engine, engine_errors = make_engine(host_python, host_error)
        secs, tl, total, speed, tm, dropped = plan_audio(engine)
        sr = engine.sr
        wav = WORK / 'narration.wav'
        track = narration_track(tl, total, sr, tm['lead'])
        write_wav(wav, track, sr)
        pre = measure_loudness(ff, wav)
        gain_db = 0.0
        if pre.get('input_i') is not None and pre['input_i'] > -70:
            # +3.01 dB once the mono voice sits in both stereo channels
            gain_db = round(min(12.0, max(0.0, TARGET_LUFS - 3.01 - pre['input_i'])), 2)
        log(f'narration loudness (mono): {pre}; gain {gain_db} dB before the -1.5 dBFS limiter')

        video = OUT / VIDEO_NAME
        pngs = [OUT / s['png'] for s in slides_info]
        log(f'encoding {len(pngs)} slides, {total:.2f}s ...')
        encode(ff, pngs, [s['dur'] for s in tl], wav, total, video, gain_db)
        info = probe(ff, fp, video)
        fmt = info.get('format', {})
        dur = float(fmt.get('duration') or 0.0)
        size = video.stat().st_size
        vs = next((s for s in info.get('streams', []) if s.get('codec_type') == 'video'), {})
        aus = next((s for s in info.get('streams', []) if s.get('codec_type') == 'audio'), {})
        log(f'video: {dur:.2f}s, {size / 1e6:.2f} MB, {vs.get("width")}x{vs.get("height")}')
        post = measure_loudness(ff, video)
        log(f'video loudness: {post}')

        thumbs = []
        for n, (idx, label) in enumerate([(0, 'title'), (2, 'slide2'), (4, 'slide4'),
                                          (7, 'closing')], 1):
            sl = tl[idx]
            name = f'thumb_{n}_{label}.png'
            grab_frame(ff, video, sl['start'] + sl['dur'] / 2, OUT / name)
            thumbs.append({'file': name, 'slide_index': idx,
                           'time_s': round(sl['start'] + sl['dur'] / 2, 2)})
        contact_sheet(ff, video, tl, OUT / 'contact_sheet.png')
        n_cues = write_captions(tl, secs, tm['lead'], OUT / 'captions.srt')

        with open(OUT / 'narration.txt', 'w', encoding='utf-8') as f:
            for sl, sents in zip(tl[1:], secs):
                st = sl['start'] + tm['lead']
                f.write(f'Slide {sl["index"]}: {SECTION_TITLES[sl["index"] - 1]} '
                        f'[{mmss(st)} - {mmss(st + sl["a_len"])}]\n')
                f.write(' '.join(s['text'] for s in sents) + '\n\n')
        with open(OUT / 'phonemes.txt', 'w', encoding='utf-8') as f:
            f.write(f'TTS: {engine.name}, voice {engine.voice}, speed {speed}, {engine.versions}\n\n')
            for k, sents in enumerate(secs, 1):
                for s in sents:
                    f.write(f'[{k}] {s["text"]}\n    input:    {s["tts_input"]}\n'
                            f'    phonemes: {s["phonemes"]}\n\n')

        # loudness of each narrated slide, from the track itself
        slide_rows, words_total, speech_s = [], 0, 0.0
        for sl in tl:
            row = dict(index=sl['index'], key=SLIDES[sl['index']][0], title=SLIDES[sl['index']][1],
                       start_s=round(sl['start'], 3), end_s=round(sl['start'] + sl['dur'], 3),
                       duration_s=round(sl['dur'], 3))
            if sl['index'] > 0:
                text = NARRATION[sl['index'] - 1]
                a0 = sl['start'] + tm['lead']
                seg = track[int(a0 * sr):int((a0 + sl['a_len']) * sr)]
                row.update(narration_start_s=round(a0, 3), narration_end_s=round(a0 + sl['a_len'], 3),
                           words=count_words(text),
                           narration_rms_dbfs_before_gain=round(
                               float(20 * np.log10(np.sqrt(np.mean(seg ** 2)) + 1e-9)), 1))
                words_total += row['words']
                speech_s += sl['a_len']
            else:
                row.update(narration_start_s=None, narration_end_s=None, words=0)
            row.update({k: v for k, v in slides_info[sl['index']].items() if k in ('png', 'font_scale', 'overflow')})
            slide_rows.append(row)

        summary.update({
            'status': 'ok',
            'video_file': VIDEO_NAME,
            'duration_s': round(dur, 3),
            'duration_mmss': mmss(dur),
            'planned_duration_s': round(total, 3),
            'under_3_min': dur < 180.0,
            'in_target_2m40_2m55': 160.0 <= dur <= 175.0,
            'file_size_bytes': size,
            'file_size_mb': round(size / 1e6, 2),
            'video': {'codec': vs.get('codec_name'), 'profile': vs.get('profile'),
                      'width': vs.get('width'), 'height': vs.get('height'),
                      'pix_fmt': vs.get('pix_fmt'), 'frame_rate': vs.get('r_frame_rate')},
            'audio': {'codec': aus.get('codec_name'), 'sample_rate': aus.get('sample_rate'),
                      'channels': aus.get('channels')},
            'tts_engine': engine.name,
            'tts_versions': engine.versions,
            'voice': engine.voice,
            'voice_language': 'American English' if engine.name.startswith('Kokoro') else 'en_US',
            'speech_speed': speed,
            'tts_errors_before_fallback': engine_errors,
            'pronunciation_hints': sorted(KOKORO_LEXICON if engine.name.startswith('Kokoro')
                                          else PIPER_RESPELL),
            'timing_settings_s': tm,
            'words_spoken': words_total,
            'speech_seconds': round(speech_s, 2),
            'words_per_minute_while_speaking': round(words_total / (speech_s / 60.0), 1) if speech_s else None,
            'slides': slide_rows,
            'loudness': dict(volume_stats(ff, video),
                             integrated_lufs=post.get('input_i'), true_peak_dbtp=post.get('input_tp'),
                             loudness_range_lu=post.get('input_lra'), target_lufs=TARGET_LUFS,
                             gain_applied_db=gain_db, narration_mono_lufs_before=pre.get('input_i')),
            'thumbnails': thumbs,
            'contact_sheet': 'contact_sheet.png',
            'captions': {'file': 'captions.srt', 'cues': n_cues},
            'repo_url_on_closing_slide': REPO_URL or None,
            'dropped_words': dropped,
            'runtime_s': round(time.time() - T0, 1),
        })
        if dropped:
            raise RuntimeError(f'the TTS dropped words it could not pronounce: {dropped}')
        if not summary['under_3_min']:
            raise RuntimeError(f'video is {dur:.2f}s, not under 3:00')
    except Exception as e:
        summary.update(status='error', error=f'{type(e).__name__}: {e}',
                       traceback=traceback.format_exc(), log_tail=LOG[-60:])
        raise
    finally:
        with open(OUT / 'summary.json', 'w', encoding='utf-8') as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)
        log(f'summary.json written (status {summary["status"]})')


if __name__ == '__main__':
    main()
