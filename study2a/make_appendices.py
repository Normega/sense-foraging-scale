"""Appendices A and B of study2a/preregistration_draft.md, from the definitions built in
radlab (scripts/sense_foraging/study2a/). Rerun if the wording changes, so the registration
quotes exactly what participants see.

Usage:  python study2a/make_appendices.py <radlab>/scripts/sense_foraging/study2a
"""
import json
import sys
from pathlib import Path

src = Path(sys.argv[1])
data = json.loads((src / 'items.json').read_text(encoding='utf-8'))
scale = json.loads((src / 'sf2a-scale.json').read_text(encoding='utf-8'))
text = {c['id']: c['question'] for p in scale['pages'] for c in p['components'] if c['type'] != 'information'}
pairs = {}
for a, b in scale['interleave']['apart_pages']:
    pairs.setdefault(a, []).append(b)
    pairs.setdefault(b, []).append(a)
src_id = {i['id']: i['source'] for i in scale['items']}

out = ['## Appendix A. The two modes', '', '**Introduction** (shown to everyone before the vignettes):', '']
out += [f'> {para}\n>' for para in data['intro']][:-1] + [f"> {data['intro'][-1]}", '']
out += [f"**Vignettes.** \"{data['stem']}\" 1 {data['anchors'][0]} to 6 {data['anchors'][1]}.", '',
        '| Id | Kind | Vignette |', '|---|---|---|']
out += [f"| {v['id']} (`vig_{v['id'][1:]}`) | {v['kind']} | {v['text']} |" for v in data['vignettes']]
out += ['', f"**Self-rating** (`vig_self`): \"{data['selfRating']}\"", '',
        '## Appendix B. The statements', '',
        'Ids in data order. R = reverse-keyed. "Study 1" = the Study 1 id of a statement still in its Study 1 wording. '
        '"Apart from" = its look-alike partners, never on the same page.', '']
for name, f in scale['factors'].items():
    out += [f'### {name}', '', f"*Definition (not shown to participants):* {f['definition']}", '',
            '| Id | Review id | Keying | Study 1 | Statement | Apart from |', '|---|---|---|---|---|---|']
    meta = {i['id']: i for i in scale['items']}
    for i in f['item_ids']:
        m = meta[i]
        out.append(f"| `{i}` | {m['source']} | {'R' if m['reverse_keyed'] else ''} | {m.get('study1_item', '')} | {text[i]} | "
                   f"{', '.join(src_id[x] for x in pairs.get(i, []))} |")
    out.append('')
out += ['### Attention checks', '', '| Id | Statement | Correct |', '|---|---|---|']
out += [f"| `{k}` | {text[k]} | {v} |" for k, v in scale['attention_checks'].items()]
out += ['', '### At the end, in Study 1 wording', '', '| Id | Statement |', '|---|---|']
out += [f"| `{e['id']}` | {e['text']} |" for e in data['end']]

p = Path(__file__).parent / 'preregistration_draft.md'
t = p.read_text(encoding='utf-8')
mark = '<!-- APPENDICES -->'
t = t[:t.index(mark) + len(mark)] + '\n\n' + '\n'.join(out) + '\n'
p.write_bytes(t.encode('utf-8'))
print(len(out), 'lines')
