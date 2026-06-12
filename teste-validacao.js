'use strict';
/* Validação numérica — replica exatamente o núcleo físico do artefato:
   mesmas constantes, dados, inicialização e integrador (leapfrog KDK, dt = 60 s).
   Verifica: (1) período orbital medido vs. Kepler, (2) deriva de energia total,
   (3) estabilidade do semieixo maior, ao longo de 600 dias simulados. */

const TAU = Math.PI * 2, DEG = Math.PI / 180;
const G = 6.6743e-20; // km³ kg⁻¹ s⁻²
const BODIES = [
  { name: 'Saturno',  mass: 5.6834e26, radius: 58232 },
  { name: 'Mimas',    mass: 3.7493e19, a: 185539,   e: 0.0196, i: 1.574,  O: 120, w: 60,  nu: 20  },
  { name: 'Encélado', mass: 1.0802e20, a: 237948,   e: 0.0047, i: 0.009,  O: 40,  w: 120, nu: 210 },
  { name: 'Tétis',    mass: 6.1745e20, a: 294619,   e: 0.0001, i: 1.091,  O: 200, w: 0,   nu: 300 },
  { name: 'Dione',    mass: 1.0955e21, a: 377396,   e: 0.0022, i: 0.028,  O: 80,  w: 200, nu: 95  },
  { name: 'Reia',     mass: 2.3065e21, a: 527108,   e: 0.0010, i: 0.331,  O: 300, w: 90,  nu: 160 },
  { name: 'Titã',     mass: 1.3452e23, a: 1221870,  e: 0.0288, i: 0.348,  O: 25,  w: 185, nu: 330 },
  { name: 'Hipérion', mass: 5.62e18,   a: 1481010,  e: 0.1230, i: 0.430,  O: 150, w: 280, nu: 70  },
  { name: 'Jápeto',   mass: 1.8056e21, a: 3560820,  e: 0.0286, i: 15.47,  O: 260, w: 50,  nu: 200 },
  { name: 'Febe',     mass: 8.292e18,  a: 12952000, e: 0.1635, i: 173.0,  O: 0,   w: 340, nu: 45  },
];
const N = BODIES.length;
const P = {
  x: new Float64Array(N), y: new Float64Array(N), z: new Float64Array(N),
  vx: new Float64Array(N), vy: new Float64Array(N), vz: new Float64Array(N),
  ax: new Float64Array(N), ay: new Float64Array(N), az: new Float64Array(N),
  m: new Float64Array(N)
};

function elementsToState(mu, a, e, i, O, w, nu) {
  i *= DEG; O *= DEG; w *= DEG; nu *= DEG;
  const p = a * (1 - e * e), r = p / (1 + e * Math.cos(nu)), h = Math.sqrt(mu * p);
  const xp = r * Math.cos(nu), yp = r * Math.sin(nu);
  const vxp = -(mu / h) * Math.sin(nu), vyp = (mu / h) * (e + Math.cos(nu));
  const cO = Math.cos(O), sO = Math.sin(O), ci = Math.cos(i), si = Math.sin(i), cw = Math.cos(w), sw = Math.sin(w);
  const R11 = cO * cw - sO * sw * ci, R12 = -cO * sw - sO * cw * ci;
  const R21 = sO * cw + cO * sw * ci, R22 = -sO * sw + cO * cw * ci;
  const R31 = sw * si, R32 = cw * si;
  const X = R11 * xp + R12 * yp, Y = R21 * xp + R22 * yp, Z = R31 * xp + R32 * yp;
  const VX = R11 * vxp + R12 * vyp, VY = R21 * vxp + R22 * vyp, VZ = R31 * vxp + R32 * vyp;
  return { x: X, y: Z, z: -Y, vx: VX, vy: VZ, vz: -VY };
}

function computeAccel() {
  P.ax.fill(0); P.ay.fill(0); P.az.fill(0);
  for (let i = 0; i < N; i++) {
    const xi = P.x[i], yi = P.y[i], zi = P.z[i];
    for (let j = i + 1; j < N; j++) {
      const dx = P.x[j] - xi, dy = P.y[j] - yi, dz = P.z[j] - zi;
      const r2 = dx * dx + dy * dy + dz * dz, inv = G / (r2 * Math.sqrt(r2));
      const fi = inv * P.m[j], fj = inv * P.m[i];
      P.ax[i] += dx * fi; P.ay[i] += dy * fi; P.az[i] += dz * fi;
      P.ax[j] -= dx * fj; P.ay[j] -= dy * fj; P.az[j] -= dz * fj;
    }
  }
}
function stepLeapfrog(dt) {
  const h = dt * 0.5;
  for (let k = 0; k < N; k++) {
    P.vx[k] += P.ax[k] * h; P.vy[k] += P.ay[k] * h; P.vz[k] += P.az[k] * h;
    P.x[k] += P.vx[k] * dt; P.y[k] += P.vy[k] * dt; P.z[k] += P.vz[k] * dt;
  }
  computeAccel();
  for (let k = 0; k < N; k++) { P.vx[k] += P.ax[k] * h; P.vy[k] += P.ay[k] * h; P.vz[k] += P.az[k] * h; }
}
function energy() {
  let E = 0;
  for (let k = 0; k < N; k++) E += 0.5 * P.m[k] * (P.vx[k] ** 2 + P.vy[k] ** 2 + P.vz[k] ** 2);
  for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
    const dx = P.x[j] - P.x[i], dy = P.y[j] - P.y[i], dz = P.z[j] - P.z[i];
    E -= G * P.m[i] * P.m[j] / Math.sqrt(dx * dx + dy * dy + dz * dz);
  }
  return E;
}
function initBodies() {
  for (let k = 0; k < N; k++) P.m[k] = BODIES[k].mass;
  P.x[0] = P.y[0] = P.z[0] = P.vx[0] = P.vy[0] = P.vz[0] = 0;
  for (let k = 1; k < N; k++) {
    const b = BODIES[k], mu = G * (BODIES[0].mass + b.mass);
    const s = elementsToState(mu, b.a, b.e, b.i, b.O, b.w, b.nu);
    P.x[k] = s.x; P.y[k] = s.y; P.z[k] = s.z; P.vx[k] = s.vx; P.vy[k] = s.vy; P.vz[k] = s.vz;
    b.kepler = TAU * Math.sqrt(b.a ** 3 / mu);
  }
  let M = 0, cx = 0, cy = 0, cz = 0, px = 0, py = 0, pz = 0;
  for (let k = 0; k < N; k++) {
    M += P.m[k];
    cx += P.m[k] * P.x[k]; cy += P.m[k] * P.y[k]; cz += P.m[k] * P.z[k];
    px += P.m[k] * P.vx[k]; py += P.m[k] * P.vy[k]; pz += P.m[k] * P.vz[k];
  }
  for (let k = 0; k < N; k++) {
    P.x[k] -= cx / M; P.y[k] -= cy / M; P.z[k] -= cz / M;
    P.vx[k] -= px / M; P.vy[k] -= py / M; P.vz[k] -= pz / M;
  }
  computeAccel();
}

/* ---------- medição ---------- */
initBodies();
const E0 = energy();
const dt = 60, DAYS = 600, STEPS = Math.round(DAYS * 86400 / dt);

// base ortonormal do plano orbital inicial de cada lua (para ângulo desenrolado)
const track = [];
for (let k = 1; k < N; k++) {
  const rx = P.x[k] - P.x[0], ry = P.y[k] - P.y[0], rz = P.z[k] - P.z[0];
  const vx = P.vx[k] - P.vx[0], vy = P.vy[k] - P.vy[0], vz = P.vz[k] - P.vz[0];
  const hx = ry * vz - rz * vy, hy = rz * vx - rx * vz, hz = rx * vy - ry * vx;
  const hn = Math.hypot(hx, hy, hz), rn = Math.hypot(rx, ry, rz);
  const e1 = [rx / rn, ry / rn, rz / rn], n = [hx / hn, hy / hn, hz / hn];
  const e2 = [n[1] * e1[2] - n[2] * e1[1], n[2] * e1[0] - n[0] * e1[2], n[0] * e1[1] - n[1] * e1[0]];
  track.push({ k, e1, e2, prev: 0, acc: 0, period: null, aMin: Infinity, aMax: -Infinity });
}
const wrapPi = a => { while (a > Math.PI) a -= TAU; while (a < -Math.PI) a += TAU; return a; };

let maxDrift = 0;
for (let s = 1; s <= STEPS; s++) {
  stepLeapfrog(dt);
  for (const t of track) {
    const k = t.k;
    const rx = P.x[k] - P.x[0], ry = P.y[k] - P.y[0], rz = P.z[k] - P.z[0];
    const ang = Math.atan2(rx * t.e2[0] + ry * t.e2[1] + rz * t.e2[2],
                           rx * t.e1[0] + ry * t.e1[1] + rz * t.e1[2]);
    const d = wrapPi(ang - t.prev); t.prev = ang;
    const before = Math.abs(t.acc); t.acc += d;
    if (t.period === null && Math.abs(t.acc) >= TAU) {
      const over = Math.abs(t.acc) - TAU;
      t.period = s * dt - dt * over / Math.abs(d || 1e-12);
    }
    if (s % 50 === 0) { // semieixo via vis-viva
      const v2 = (P.vx[k] - P.vx[0]) ** 2 + (P.vy[k] - P.vy[0]) ** 2 + (P.vz[k] - P.vz[0]) ** 2;
      const r = Math.hypot(rx, ry, rz), mu = G * (BODIES[0].mass + BODIES[k].mass);
      const a = 1 / (2 / r - v2 / mu);
      if (a < t.aMin) t.aMin = a; if (a > t.aMax) t.aMax = a;
    }
  }
  if (s % 720 === 0) { // energia a cada 12 h simuladas
    const drift = Math.abs((energy() - E0) / E0);
    if (drift > maxDrift) maxDrift = drift;
  }
}

/* ---------- relatório ---------- */
const pad = (s, n) => String(s).padEnd(n);
console.log('VALIDAÇÃO — ' + DAYS + ' dias simulados, dt = ' + dt + ' s, ' + STEPS.toLocaleString('pt-BR') + ' passos, N = ' + N + ' corpos\n');
console.log(pad('Lua', 10) + pad('P Kepler (d)', 14) + pad('P medido (d)', 14) + pad('erro %', 10) + pad('Δa %', 10) + 'status');
let ok = true;
for (const t of track) {
  const b = BODIES[t.k];
  const pk = b.kepler / 86400, pm = t.period / 86400;
  const err = Math.abs(pm - pk) / pk * 100;
  const da = (t.aMax - t.aMin) / b.a * 100;
  const lim = b.name === 'Hipérion' ? 1.5 : 0.5; // ressonância 4:3 com Titã perturba mais
  const pass = err < lim && da < (b.name === 'Hipérion' ? 3 : 1);
  if (!pass) ok = false;
  console.log(pad(b.name, 10) + pad(pk.toFixed(3), 14) + pad(pm.toFixed(3), 14) +
    pad(err.toFixed(4), 10) + pad(da.toFixed(4), 10) + (pass ? 'OK' : 'FALHOU'));
}
console.log('\nDeriva máxima de energia total: ' + maxDrift.toExponential(2) +
  ' (relativa)  →  ' + (maxDrift < 1e-5 ? 'OK' : 'FALHOU'));
if (maxDrift >= 1e-5) ok = false;
console.log(ok ? '\nRESULTADO: TODOS OS TESTES PASSARAM' : '\nRESULTADO: HÁ FALHAS');
process.exit(ok ? 0 : 1);
