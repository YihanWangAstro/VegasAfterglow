"use client";

import dynamic from "next/dynamic";

const ClientPage = dynamic(() => import("./client-page"), {
  ssr: false,
  loading: () => (
    <main className="page">
      <section className="panel">
        <h1>VegasAfterglow Webtool</h1>
        <p className="muted">Loading client UI...</p>
      </section>
    </main>
  ),
});

export default function Page() {
  return <ClientPage />;
}
