export function HelpHint({ text, ariaLabel }: { text: string; ariaLabel: string }) {
  return (
    <span className="sb-help" tabIndex={0} role="button" aria-label={ariaLabel}>
      ?
      <span className="sb-help-tip">{text}</span>
    </span>
  );
}
