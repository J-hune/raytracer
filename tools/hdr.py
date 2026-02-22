"""Lecture et ecriture de fichiers Radiance .hdr (RGBE).

L'ecriture produit des scanlines RLE. Ecrire les donnees a plat parait plus
simple mais stb_image, qu'utilise le moteur, bascule sur son chemin RLE des que
la largeur est dans [8, 32768) et son repli non-RLE ne remplit qu'une ligne sur
deux, le reste sortant a zero.
"""
import numpy as np


def read(path):
    """Renvoie un tableau (H, W, 3) de flottants lineaires."""
    with open(path, "rb") as handle:
        while handle.readline().strip():
            pass
        dims = handle.readline().split()
        height, width = int(dims[1]), int(dims[3])

        raw = np.zeros((height, width, 4), np.uint8)
        for y in range(height):
            head = handle.read(4)
            if head[0] == 2 and head[1] == 2:
                for channel in range(4):
                    x = 0
                    while x < width:
                        n = handle.read(1)[0]
                        if n > 128:
                            raw[y, x:x + n - 128, channel] = handle.read(1)[0]
                            x += n - 128
                        else:
                            raw[y, x:x + n, channel] = np.frombuffer(
                                handle.read(n), np.uint8)
                            x += n
            else:
                rest = np.frombuffer(handle.read(4 * width - 4), np.uint8)
                raw[y, 0] = np.frombuffer(head, np.uint8)
                raw[y, 1:] = rest.reshape(width - 1, 4)

    exponent = raw[:, :, 3].astype(np.int32)
    scale = np.where(exponent > 0, np.exp2(exponent - 136.0), 0.0)
    return raw[:, :, :3].astype(np.float64) * scale[:, :, None]


def _encode_scanline(row):
    """RLE Radiance sur un canal d'une scanline."""
    out = bytearray()
    width = len(row)
    i = 0
    while i < width:
        run = 1
        while (i + run < width and run < 127 and row[i + run] == row[i]):
            run += 1
        if run >= 4:
            out.append(128 + run)
            out.append(int(row[i]))
            i += run
            continue
        start = i
        while i < width and i - start < 128:
            ahead = 1
            while (i + ahead < width and ahead < 4 and row[i + ahead] == row[i]):
                ahead += 1
            if ahead >= 4:                      # une repetition commence: on coupe
                break
            i += 1
        out.append(i - start)
        out.extend(row[start:i].tobytes())
    return out


def write(path, rgb, dither=True, seed=0):
    """Ecrit un tableau (H, W, 3) de flottants lineaires."""
    height, width = rgb.shape[:2]
    peak = rgb.max(axis=2)
    raw = np.zeros((height, width, 4), np.uint8)
    live = peak > 1e-32
    mantissa, exponent = np.frexp(peak[live])
    raw[live, 3] = (exponent + 128).astype(np.uint8)
    scaled = rgb[live] * (mantissa * 256.0 / peak[live])[:, None]
    if dither:
        # La mantisse n'a que 8 bits: sans tramage un degrade doux quantifie en
        # paliers horizontaux visibles.
        scaled = scaled + np.random.default_rng(seed).random(scaled.shape)
    else:
        scaled = scaled + 0.5
    raw[live, 0:3] = np.clip(np.floor(scaled), 0, 255).astype(np.uint8)

    with open(path, "wb") as handle:
        handle.write(b"#?RADIANCE\nFORMAT=32-bit_rle_rgbe\n\n")
        handle.write(b"-Y %d +X %d\n" % (height, width))
        for y in range(height):
            handle.write(bytes([2, 2, (width >> 8) & 0xFF, width & 0xFF]))
            for channel in range(4):
                handle.write(_encode_scanline(raw[y, :, channel]))
