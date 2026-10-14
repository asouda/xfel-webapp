(() => {
  const canvas = document.getElementById("field");
  if (!canvas) return;
  const ctx = canvas.getContext("2d");
  const reduce = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  const colors = [
    [90, 150, 255],
    [70, 120, 255],
    [150, 90, 255],
    [190, 80, 255],
    [120, 160, 255],
  ];
  let particles = [];

  function resize() {
    canvas.width = window.innerWidth;
    canvas.height = window.innerHeight;
    seed();
    draw();
  }

  function bell() {
    return (Math.random() + Math.random() + Math.random()) / 3;
  }

  function seed() {
    const count = Math.round(Math.min(110, Math.max(50, canvas.width * canvas.height / 18000)));
    particles = [];
    for (let i = 0; i < count; i += 1) {
      const color = colors[i % colors.length];
      const large = Math.random() < 0.08;
      particles.push({
        x: bell() * canvas.width,
        y: bell() * canvas.height,
        r: large ? 18 + Math.random() * 46 : 1.2 + Math.random() * 3.4,
        vx: (Math.random() - 0.5) * (large ? 0.08 : 0.22),
        vy: (Math.random() - 0.5) * (large ? 0.08 : 0.22),
        a: large ? 0.05 + Math.random() * 0.08 : 0.25 + Math.random() * 0.55,
        color,
      });
    }
  }

  function draw() {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    for (const p of particles) {
      ctx.beginPath();
      ctx.fillStyle = `rgba(${p.color[0]}, ${p.color[1]}, ${p.color[2]}, ${p.a})`;
      ctx.shadowColor = `rgba(${p.color[0]}, ${p.color[1]}, ${p.color[2]}, ${Math.min(0.9, p.a + 0.2)})`;
      ctx.shadowBlur = p.r > 10 ? 40 : 16;
      ctx.arc(p.x, p.y, p.r, 0, Math.PI * 2);
      ctx.fill();
    }
    ctx.shadowBlur = 0;
  }

  function tick() {
    for (const p of particles) {
      p.x += p.vx;
      p.y += p.vy;
      if (p.x < -40) p.x = canvas.width + 40;
      if (p.x > canvas.width + 40) p.x = -40;
      if (p.y < -40) p.y = canvas.height + 40;
      if (p.y > canvas.height + 40) p.y = -40;
    }
    draw();
    requestAnimationFrame(tick);
  }

  window.addEventListener("resize", resize);
  resize();
  if (!reduce) requestAnimationFrame(tick);
})();
