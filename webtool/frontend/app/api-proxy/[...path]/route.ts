import { NextRequest, NextResponse } from "next/server";

export const runtime = "nodejs";
export const dynamic = "force-dynamic";

function stripTrailingSlash(url: string): string {
  return url.replace(/\/+$/, "");
}

function backendCandidates(): string[] {
  const configured = process.env.BACKEND_URL || process.env.NEXT_PUBLIC_API_URL || "http://127.0.0.1:8000";
  const primary = stripTrailingSlash(configured);
  const list = [primary];
  if (primary.includes("localhost")) {
    list.push(primary.replace("localhost", "127.0.0.1"));
  }
  if (primary.includes("127.0.0.1")) {
    list.push(primary.replace("127.0.0.1", "localhost"));
  }
  return Array.from(new Set(list));
}

type RouteContext = { params: Promise<{ path: string[] }> };

async function forward(request: NextRequest, context: RouteContext) {
  try {
    const resolved = await context.params;
    const pathParts = Array.isArray(resolved?.path)
      ? resolved.path
      : typeof resolved?.path === "string"
        ? [resolved.path]
        : [];
    const path = pathParts.join("/");
    const query = request.nextUrl.search;

    const method = request.method;
    const includeBody = !["GET", "HEAD"].includes(method);
    const rawBody = includeBody ? await request.text() : "";
    const body = rawBody.length > 0 ? rawBody : undefined;

    let lastError: unknown = null;
    const tried = backendCandidates();

    for (const base of tried) {
      const target = `${base}/api/${path}${query}`;
      try {
        const headers: Record<string, string> = {};
        const contentType = request.headers.get("content-type");
        if (contentType) {
          headers["content-type"] = contentType;
        }

        const response = await fetch(target, {
          method,
          headers,
          body,
          redirect: "manual",
        });

        const responseBody = await response.arrayBuffer();
        const responseType = response.headers.get("content-type") ?? "application/json";
        return new NextResponse(responseBody, {
          status: response.status,
          statusText: response.statusText,
          headers: { "content-type": responseType },
        });
      } catch (err) {
        lastError = err;
        // eslint-disable-next-line no-console
        console.error(`API proxy failed for ${target}:`, err);
      }
    }

    const message = lastError instanceof Error ? lastError.message : "Unknown proxy error";
    return NextResponse.json(
      {
        error: "Backend unreachable",
        detail: message,
        tried,
      },
      { status: 502 },
    );
  } catch (err) {
    const message = err instanceof Error ? err.message : "Unknown proxy crash";
    const stack = err instanceof Error ? err.stack : undefined;
    return NextResponse.json(
      {
        error: "Proxy handler crashed",
        detail: message,
        stack,
      },
      { status: 500 },
    );
  }
}

export async function GET(request: NextRequest, context: RouteContext) {
  return forward(request, context);
}

export async function POST(request: NextRequest, context: RouteContext) {
  return forward(request, context);
}

export async function PUT(request: NextRequest, context: RouteContext) {
  return forward(request, context);
}

export async function PATCH(request: NextRequest, context: RouteContext) {
  return forward(request, context);
}

export async function DELETE(request: NextRequest, context: RouteContext) {
  return forward(request, context);
}

export async function OPTIONS(request: NextRequest, context: RouteContext) {
  return forward(request, context);
}
