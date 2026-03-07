#!/usr/bin/env node
import assert from "node:assert/strict";
import fs from "node:fs";
import { Buffer } from "node:buffer";

function luminosityDistanceMpcFromRedshift(z) {
  const H0 = 67.4;
  const C = 299792.458;
  const zSafe = Math.max(0, z);
  return (C / H0) * zSafe * (1 + zSafe);
}

function redshiftFromLuminosityDistanceMpc(dL) {
  const H0 = 67.4;
  const C = 299792.458;
  const dSafe = Math.max(0, dL);
  const x = (dSafe * H0) / C;
  return (-1 + Math.sqrt(1 + 4 * x)) / 2;
}

function encodeState(snapshot) {
  const payload = JSON.stringify({ v: 1, snapshot });
  const bytes = new TextEncoder().encode(payload);
  let binary = "";
  for (let i = 0; i < bytes.length; i += 1) binary += String.fromCharCode(bytes[i]);
  return Buffer.from(binary, "binary")
    .toString("base64")
    .replace(/\+/g, "-")
    .replace(/\//g, "_")
    .replace(/=+$/g, "");
}

function decodeState(raw) {
  const normalized = raw.replace(/-/g, "+").replace(/_/g, "/");
  const padded = normalized + "=".repeat((4 - (normalized.length % 4)) % 4);
  const binary = Buffer.from(padded, "base64").toString("binary");
  const bytes = Uint8Array.from(binary, (c) => c.charCodeAt(0));
  return JSON.parse(new TextDecoder().decode(bytes));
}

const NO_DECODE = Symbol("no-decode");
const DTYPE_BYTES = {
  f4: 4,
  f8: 8,
  i1: 1,
  u1: 1,
  i2: 2,
  u2: 2,
  i4: 4,
  u4: 4,
};

function parseShape(raw) {
  if (Array.isArray(raw)) {
    return raw
      .map((value) => Number(value))
      .filter((value) => Number.isInteger(value) && value > 0);
  }
  if (typeof raw === "string") {
    return raw
      .split(",")
      .map((part) => Number(part.trim()))
      .filter((value) => Number.isInteger(value) && value > 0);
  }
  return [];
}

function reshapeFlat(values, shape) {
  if (shape.length <= 1) return values;
  const expected = shape.reduce((acc, dim) => acc * dim, 1);
  if (expected !== values.length) return values;

  const build = (offset, dims) => {
    if (dims.length === 1) {
      const size = dims[0];
      return [values.slice(offset, offset + size), offset + size];
    }
    const [head, ...rest] = dims;
    const out = new Array(head);
    let cursor = offset;
    for (let i = 0; i < head; i += 1) {
      const [chunk, next] = build(cursor, rest);
      out[i] = chunk;
      cursor = next;
    }
    return [out, cursor];
  };

  return build(0, shape)[0];
}

function decodeBinaryNode(value) {
  if (!value || typeof value !== "object") return NO_DECODE;
  const dtype = typeof value.dtype === "string" ? value.dtype : null;
  const bdata = typeof value.bdata === "string" ? value.bdata : null;
  if (!dtype || !bdata || !(dtype in DTYPE_BYTES)) return NO_DECODE;

  let bytes;
  try {
    bytes = Buffer.from(bdata, "base64");
  } catch {
    return NO_DECODE;
  }

  const byteWidth = DTYPE_BYTES[dtype];
  if (bytes.byteLength === 0 || bytes.byteLength % byteWidth !== 0) return NO_DECODE;

  try {
    const buffer = bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength);
    let flat;
    switch (dtype) {
      case "f4":
        flat = Array.from(new Float32Array(buffer), Number);
        break;
      case "f8":
        flat = Array.from(new Float64Array(buffer), Number);
        break;
      case "i1":
        flat = Array.from(new Int8Array(buffer), Number);
        break;
      case "u1":
        flat = Array.from(new Uint8Array(buffer), Number);
        break;
      case "i2":
        flat = Array.from(new Int16Array(buffer), Number);
        break;
      case "u2":
        flat = Array.from(new Uint16Array(buffer), Number);
        break;
      case "i4":
        flat = Array.from(new Int32Array(buffer), Number);
        break;
      case "u4":
        flat = Array.from(new Uint32Array(buffer), Number);
        break;
      default:
        return NO_DECODE;
    }
    return reshapeFlat(flat, parseShape(value.shape));
  } catch {
    return NO_DECODE;
  }
}

function decodePlotlyFigureBinaryInPlace(figure) {
  if (!figure || typeof figure !== "object") return;
  const stack = [figure];

  while (stack.length > 0) {
    const current = stack.pop();
    if (!current || typeof current !== "object") continue;

    if (Array.isArray(current)) {
      for (let i = 0; i < current.length; i += 1) {
        const item = current[i];
        if (!item || typeof item !== "object") continue;
        const decoded = decodeBinaryNode(item);
        if (decoded !== NO_DECODE) {
          current[i] = decoded;
        } else {
          stack.push(item);
        }
      }
      continue;
    }

    for (const key of Object.keys(current)) {
      const item = current[key];
      if (!item || typeof item !== "object") continue;
      const decoded = decodeBinaryNode(item);
      if (decoded !== NO_DECODE) {
        current[key] = decoded;
      } else {
        stack.push(item);
      }
    }
  }
}

function makeF4Binary(values, shape) {
  const typed = Float32Array.from(values);
  const bytes = Buffer.from(typed.buffer.slice(typed.byteOffset, typed.byteOffset + typed.byteLength));
  return {
    dtype: "f4",
    bdata: bytes.toString("base64"),
    shape,
  };
}

function assertArrayClose(actual, expected, eps = 1e-6) {
  assert.equal(actual.length, expected.length, "array length mismatch");
  for (let i = 0; i < actual.length; i += 1) {
    assert.ok(Math.abs(actual[i] - expected[i]) < eps, `array mismatch at ${i}`);
  }
}

function main() {
  const z = 0.022;
  const dL = luminosityDistanceMpcFromRedshift(z);
  const zBack = redshiftFromLuminosityDistanceMpc(dL);
  assert.ok(Math.abs(z - zBack) < 1e-12, "distance z<->dL conversion drift");

  const snapshot = {
    mode: "lightcurve",
    distance_linked: true,
    distance_driver: "dL",
    shared: { d_L_mpc: 100, z: 0.022 },
  };
  const encoded = encodeState(snapshot);
  const decoded = decodeState(encoded);
  assert.equal(decoded.v, 1, "url state version mismatch");
  assert.equal(decoded.snapshot.mode, snapshot.mode, "url state mode mismatch");

  const constantsText = fs.readFileSync("app/lib/constants.tsx", "utf8");
  assert.ok(constantsText.includes("kind: \"log\""), "expected schema-style slider spec not found");
  assert.ok(constantsText.includes("kind: \"linear\""), "expected schema-style slider spec not found");

  const figure1d = {
    data: [
      {
        x: makeF4Binary([1, 2, 3], [3]),
        y: makeF4Binary([4, 5, 6], "3"),
      },
    ],
  };
  decodePlotlyFigureBinaryInPlace(figure1d);
  assert.ok(Array.isArray(figure1d.data[0].x), "1D x should decode to array");
  assert.ok(Array.isArray(figure1d.data[0].y), "1D y should decode to array");
  assertArrayClose(figure1d.data[0].x, [1, 2, 3]);
  assertArrayClose(figure1d.data[0].y, [4, 5, 6]);

  const figure2d = {
    data: [
      {
        z: makeF4Binary([1, 2, 3, 4, 5, 6], "2, 3"),
      },
    ],
  };
  decodePlotlyFigureBinaryInPlace(figure2d);
  assert.deepEqual(figure2d.data[0].z, [[1, 2, 3], [4, 5, 6]], "2D z decode failed for string shape");

  const figureFrames = {
    data: [{ z: makeF4Binary([1, 2, 3, 4], "2,2") }],
    frames: [
      { data: [{ z: makeF4Binary([10, 20, 30, 40], "2, 2") }] },
      { data: [{ z: makeF4Binary([100, 200, 300, 400], [2, 2]) }] },
    ],
  };
  decodePlotlyFigureBinaryInPlace(figureFrames);
  assert.deepEqual(figureFrames.frames[0].data[0].z, [[10, 20], [30, 40]], "frame 0 z decode failed");
  assert.deepEqual(figureFrames.frames[1].data[0].z, [[100, 200], [300, 400]], "frame 1 z decode failed");

  const badNode = { dtype: "f4", bdata: "!!!not-base64!!!", shape: "2,2" };
  const figureBad = { data: [{ x: badNode }] };
  decodePlotlyFigureBinaryInPlace(figureBad);
  assert.equal(figureBad.data[0].x, badNode, "invalid binary node should remain unchanged");

  console.log("frontend smoke ok");
}

main();
