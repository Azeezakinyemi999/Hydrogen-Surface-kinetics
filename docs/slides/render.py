"""Render a .pptx to PNG images through Microsoft PowerPoint (macOS), then pdftoppm.

PowerPoint is used because Calibri and Cambria exist only inside PowerPoint.app on this
Mac: any other renderer substitutes fonts, and line breaks would not match the real deck.

    python render.py deck.pptx out_dir [--dpi 80] [--first N --last M]
"""
import argparse
import os
import shutil
import subprocess
import time

CONTAINER = os.path.expanduser('~/Library/Containers/com.microsoft.Powerpoint/Data/deck_render')


def held(path):
    """True while PowerPoint has the file open (it keeps a ~$ lock file next to it)."""
    folder, base = os.path.split(path)
    return os.path.exists(os.path.join(folder, '~$' + base))


def clean_old_copies(name):
    """Remove copies left by earlier failed renders, but none that PowerPoint still holds."""
    for f in os.listdir(CONTAINER):
        p = os.path.join(CONTAINER, f)
        if f.startswith(name + '_render_') and f.endswith(('.pptx', '.pdf')) and not held(p):
            os.remove(p)


def powerpoint_running():
    return subprocess.run(['pgrep', '-f', 'MacOS/Microsoft PowerPoint'], capture_output=True).returncode == 0


def powerpoint_answers(seconds=20):
    """True if a running PowerPoint answers a harmless question within `seconds`."""
    q = subprocess.run(['osascript', '-e', f'with timeout of {seconds} seconds',
                        '-e', 'tell application "Microsoft PowerPoint" to get count of presentations',
                        '-e', 'end timeout'], capture_output=True, text=True)
    return q.returncode == 0


def render(pptx, out_dir, dpi=80, first=None, last=None):
    """Export through PowerPoint without touching any other open presentation: the deck is
    copied under a unique name, opened, exported and closed by that name only. If anything
    fails, the copy is closed by name; its file is kept while PowerPoint still holds it.
    A PowerPoint started here is quit afterwards when no other presentation is open, so it
    never sits on its start screen (which blocks scripted 'open' requests)."""
    os.makedirs(CONTAINER, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)
    name = os.path.splitext(os.path.basename(pptx))[0]
    was_running = powerpoint_running()
    if was_running and not powerpoint_answers():
        raise RuntimeError('PowerPoint is running but does not answer: a window there is probably waiting '
                           '(start screen, recovery, sign-in or a dialog). Dismiss it, or quit PowerPoint, and retry.')
    # A PowerPoint started here is quit while the copy is still open: closing the last presentation
    # first would bring up the start screen, which swallows the quit and later 'open' requests.
    finish = '''
    close pres saving no''' if was_running else '''
    if (count of presentations) is 1 then
        quit saving no
    else
        close pres saving no
    end if'''
    clean_old_copies(name)
    tag = f'{name}_render_{os.getpid()}_{int(time.time())}'
    src = os.path.join(CONTAINER, tag + '.pptx')
    pdf = os.path.join(CONTAINER, tag + '.pdf')
    shutil.copy(pptx, src)
    script = f'''
with timeout of 900 seconds
tell application "Microsoft PowerPoint"
    open POSIX file "{src}"
    set found to false
    repeat 120 times
        try
            set pres to presentation "{tag}.pptx"
            set found to true
            exit repeat
        end try
        delay 0.5
    end repeat
    if not found then error "the copy {tag}.pptx did not open within 60 s"
    try
        save pres in POSIX file "{pdf}" as save as PDF
    on error msg
        close pres saving no
        error msg
    end try
{finish}
end tell
end timeout
'''
    done = subprocess.run(['osascript', '-e', script], capture_output=True, text=True)
    if not held(src):
        os.remove(src)
    if done.returncode != 0 and os.path.exists(pdf) and os.path.getsize(pdf) > 0:
        print(f'warning: exported, but PowerPoint reported: {done.stderr.strip()}')
    elif done.returncode != 0:
        kept = f' The copy {src} is still open in PowerPoint: close it without saving.' if held(src) else ''
        raise RuntimeError(f'PowerPoint export failed: {done.stderr.strip()}.{kept}')
    for f in os.listdir(out_dir):
        if f.startswith(name + '-') and f.endswith('.png'):
            os.remove(os.path.join(out_dir, f))
    cmd = ['pdftoppm', '-png', '-r', str(dpi)]
    if first:
        cmd += ['-f', str(first)]
    if last:
        cmd += ['-l', str(last)]
    subprocess.run(cmd + [pdf, os.path.join(out_dir, name)], check=True)
    shutil.move(pdf, os.path.join(out_dir, name + '.pdf'))
    return sorted(f for f in os.listdir(out_dir) if f.endswith('.png'))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('pptx')
    ap.add_argument('out_dir')
    ap.add_argument('--dpi', type=int, default=80)
    ap.add_argument('--first', type=int)
    ap.add_argument('--last', type=int)
    a = ap.parse_args()
    for f in render(a.pptx, a.out_dir, a.dpi, a.first, a.last):
        print(f)
