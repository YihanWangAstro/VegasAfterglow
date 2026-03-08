import type { Metadata, Viewport } from "next";
import Script from "next/script";
import "./globals.css";

const SITE_URL = process.env.NEXT_PUBLIC_SITE_URL?.replace(/\/+$/, "") || "https://vegasafterglow.com";

export const metadata: Metadata = {
  title: "VegasAfterglow",
  description: "VegasAfterglow FastAPI + Next.js webtool",
  keywords: [
    "gamma-ray burst",
    "GRB afterglow",
    "relativistic jet",
    "afterglow light curve",
    "GRB spectrum",
    "afterglow sky image",
    "ISM wind medium",
    "multi-band astrophysics",
  ],
  metadataBase: new URL(SITE_URL),
  alternates: {
    canonical: "/",
  },
};

export const viewport: Viewport = {
  width: "device-width",
  initialScale: 1,
  maximumScale: 1,
  userScalable: false,
  viewportFit: "cover",
};

export default function RootLayout({
  children,
}: Readonly<{ children: React.ReactNode }>) {
  return (
    <html lang="en" suppressHydrationWarning>
      <head>
        <Script
          src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js"
          strategy="afterInteractive"
        />
      </head>
      <body suppressHydrationWarning>{children}</body>
    </html>
  );
}
