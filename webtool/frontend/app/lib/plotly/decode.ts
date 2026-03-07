import type { RunResponse } from "../types";

type NumericTypedArray =
  | Float32Array
  | Float64Array
  | Int8Array
  | Uint8Array
  | Int16Array
  | Uint16Array
  | Int32Array
  | Uint32Array;

type TypedArrayCtor = {
  readonly BYTES_PER_ELEMENT: number;
  new (buffer: ArrayBufferLike, byteOffset?: number, length?: number): NumericTypedArray;
};

type BinaryNode = {
  dtype?: unknown;
  bdata?: unknown;
  shape?: unknown;
};

const DTYPE_CTORS: Record<string, TypedArrayCtor> = {
  f4: Float32Array,
  f8: Float64Array,
  i1: Int8Array,
  u1: Uint8Array,
  i2: Int16Array,
  u2: Uint16Array,
  i4: Int32Array,
  u4: Uint32Array,
};

const NO_DECODE = Symbol("plotly-no-decode");

function base64ToBytes(base64: string): Uint8Array | null {
  try {
    if (typeof atob === "function") {
      const binary = atob(base64);
      const bytes = new Uint8Array(binary.length);
      for (let i = 0; i < binary.length; i += 1) {
        bytes[i] = binary.charCodeAt(i);
      }
      return bytes;
    }
    if (typeof Buffer !== "undefined") {
      return new Uint8Array(Buffer.from(base64, "base64"));
    }
  } catch {
    return null;
  }
  return null;
}

function parseShape(raw: unknown): number[] {
  if (Array.isArray(raw)) {
    const parsed = raw
      .map((value) => Number(value))
      .filter((value) => Number.isInteger(value) && value > 0);
    return parsed;
  }
  if (typeof raw === "string") {
    const parts = raw
      .split(",")
      .map((part) => Number(part.trim()))
      .filter((value) => Number.isInteger(value) && value > 0);
    if (parts.length > 0) return parts;
  }
  return [];
}

function reshapeFlat(values: number[], shape: number[]): unknown {
  if (shape.length <= 1) return values;
  if (shape.some((dim) => !Number.isInteger(dim) || dim <= 0)) return values;

  const expectedSize = shape.reduce((acc, dim) => acc * dim, 1);
  if (expectedSize !== values.length) return values;

  const makeSegment = (offset: number, dims: number[]): [unknown, number] => {
    if (dims.length === 1) {
      const size = dims[0];
      const next = values.slice(offset, offset + size);
      return [next, offset + size];
    }

    const [head, ...rest] = dims;
    const out = new Array(head);
    let cursor = offset;
    for (let i = 0; i < head; i += 1) {
      const [segment, nextCursor] = makeSegment(cursor, rest);
      out[i] = segment;
      cursor = nextCursor;
    }
    return [out, cursor];
  };

  return makeSegment(0, shape)[0];
}

function decodeBinaryNode(node: BinaryNode): unknown | typeof NO_DECODE {
  const dtype = typeof node.dtype === "string" ? node.dtype : null;
  const bdata = typeof node.bdata === "string" ? node.bdata : null;
  if (!dtype || !bdata) return NO_DECODE;

  const ctor = DTYPE_CTORS[dtype];
  if (!ctor) return NO_DECODE;

  const bytes = base64ToBytes(bdata);
  if (!bytes) return NO_DECODE;
  if (bytes.byteLength % ctor.BYTES_PER_ELEMENT !== 0) return NO_DECODE;

  const alignedBytes =
    bytes.byteOffset % ctor.BYTES_PER_ELEMENT === 0 ? bytes : new Uint8Array(bytes);

  try {
    const typed = new ctor(
      alignedBytes.buffer,
      alignedBytes.byteOffset,
      alignedBytes.byteLength / ctor.BYTES_PER_ELEMENT,
    );
    const flat = Array.from(typed, (value) => Number(value));
    const shape = parseShape(node.shape);
    return reshapeFlat(flat, shape);
  } catch {
    return NO_DECODE;
  }
}

export function decodePlotlyFigureBinaryInPlace(figure: RunResponse["figure"]): void {
  if (!figure || typeof figure !== "object") return;

  const stack: unknown[] = [figure];
  while (stack.length > 0) {
    const current = stack.pop();
    if (!current || typeof current !== "object") continue;

    if (Array.isArray(current)) {
      for (let i = 0; i < current.length; i += 1) {
        const value = current[i];
        if (value && typeof value === "object") {
          const decoded = decodeBinaryNode(value as BinaryNode);
          if (decoded !== NO_DECODE) {
            current[i] = decoded;
          } else {
            stack.push(value);
          }
        }
      }
      continue;
    }

    const record = current as Record<string, unknown>;
    for (const key of Object.keys(record)) {
      const value = record[key];
      if (!value || typeof value !== "object") continue;

      const decoded = decodeBinaryNode(value as BinaryNode);
      if (decoded !== NO_DECODE) {
        record[key] = decoded;
      } else {
        stack.push(value);
      }
    }
  }
}
