import assert from "node:assert/strict";
import { test } from "node:test";
import { Box2D } from "../box2d.js";

const { b2Vec2, b2Math } = Box2D.Common.Math;
const { b2Body, b2BodyDef, b2FilterData, b2FixtureDef, b2World } = Box2D.Dynamics;
const { b2CircleShape, b2PolygonShape } = Box2D.Collision.Shapes;

type World = InstanceType<typeof b2World>;
type Body = InstanceType<typeof b2Body>;
type Fixture = InstanceType<typeof Box2D.Dynamics.b2Fixture>;

function createBox(world: World, type: number, x0: number, y0: number, x1: number, y1: number): Body {
  const fixtureDef = new b2FixtureDef();
  fixtureDef.density = 1.0;
  fixtureDef.friction = 0.05;
  fixtureDef.restitution = 0.0;
  fixtureDef.shape = b2PolygonShape.AsArray([b2Vec2.Make(x0, y0), b2Vec2.Make(x1, y0), b2Vec2.Make(x1, y1), b2Vec2.Make(x0, y1)], 4);

  const bodyDef = new b2BodyDef();
  bodyDef.type = type;

  const body = world.CreateBody(bodyDef);
  assert.ok(body);
  body.CreateFixture(fixtureDef);
  return body;
}

function createBall(world: World, x: number, y: number, angle: number): Body {
  const fixtureDef = new b2FixtureDef();
  fixtureDef.density = 0.01;
  fixtureDef.friction = 0.0;
  fixtureDef.restitution = 1.0;
  fixtureDef.shape = new b2CircleShape(0.25);

  const bodyDef = new b2BodyDef();
  bodyDef.type = b2Body.b2_dynamicBody;
  bodyDef.fixedRotation = true;
  bodyDef.bullet = true;

  const body = world.CreateBody(bodyDef);
  assert.ok(body);
  body.CreateFixture(fixtureDef);
  body.SetPosition(b2Vec2.Make(x, y));
  body.SetLinearVelocity(b2Vec2.Make(b2Math.Cos(angle) * 18, b2Math.Sin(angle) * 18));
  return body;
}

function createCorner(): World {
  const world = new b2World(b2Vec2.Make(0, 0), true);
  createBox(world, b2Body.b2_staticBody, -0.4, -0.4, 20.4, 0.4);
  createBox(world, b2Body.b2_staticBody, -0.4, -0.4, 0.4, 20.4);
  return world;
}

function createScene() {
  const world = createCorner();
  const box = createBox(world, b2Body.b2_dynamicBody, 4, 4, 5.6, 5.2);
  const balls = [];
  for (let i = 0; i < 8; i++) {
    balls.push(createBall(world, 2 + i, 8, -2.0 - (i * 0.1)));
  }
  return { world, bodies: [box, ...balls] };
}

function stepScene(scene: ReturnType<typeof createScene>, step: number) {
  const box = scene.bodies[0];
  box.SetLinearVelocity(b2Vec2.Make(b2Math.Cos(box.GetAngle()) * -6, b2Math.Sin(box.GetAngle()) * -6));
  box.SetAngularVelocity(step % 40 < 20 ? 2 : -2);
  scene.world.Step((32 + ((step * 7) % 4)) / 1000, 10, 10);
}

function readState(bodies: Body[]): number[] {
  const state = [];
  for (const body of bodies) {
    state.push(body.GetPosition().x, body.GetPosition().y, body.GetAngle(), body.GetLinearVelocity().x, body.GetLinearVelocity().y);
  }
  return state;
}

test("time of impact keeps the bounce off the first wall in an inner corner", () => {
  const world = createCorner();
  const ball = createBall(world, 6, 3.864, b2Math.Atan2(0.4 - 3.864, 0.4 - 6));
  const initialVelocity = ball.GetLinearVelocity().Copy();

  // The ball ends the tenth step inside the first wall's contact band and reaches the second wall in the eleventh.
  for (const milliseconds of [34, 32, 34, 35, 35, 33, 34, 34, 33, 35, 33, 33]) {
    world.Step(milliseconds / 1000, 10, 10);
  }

  const velocity = ball.GetLinearVelocity();
  assert.ok(Math.abs(velocity.x + initialVelocity.x) < 1e-9, `x velocity ${velocity.x}`);
  assert.ok(Math.abs(velocity.y + initialVelocity.y) < 1e-9, `y velocity ${velocity.y}`);
});

test("SetFilterData removes contacts the new filter rejects", () => {
  const world = new b2World(b2Vec2.Make(0, 0), true);
  createBox(world, b2Body.b2_staticBody, 0, 0, 10, 1);
  const box = createBox(world, b2Body.b2_dynamicBody, 4, 0.9, 6, 2);
  world.Step(1 / 30, 10, 10);
  assert.notEqual(world.GetContactList(), null);

  const filter = new b2FilterData();
  filter.maskBits = 0;
  (box.GetFixtureList() as Fixture).SetFilterData(filter);
  world.Step(1 / 30, 10, 10);
  assert.equal(world.GetContactList(), null);
});

test("contact count follows destroyed contacts", () => {
  const world = new b2World(b2Vec2.Make(0, 0), true);
  createBox(world, b2Body.b2_staticBody, 0, 0, 10, 1);
  const box = createBox(world, b2Body.b2_dynamicBody, 4, 0.9, 6, 2);
  world.Step(1 / 30, 10, 10);
  assert.equal(world.GetContactCount(), 1);

  world.DestroyBody(box);
  assert.equal(world.GetContactCount(), 0);
});

test("worlds stepped side by side match a world stepped alone", () => {
  const alone = createScene();
  for (let step = 0; step < 120; step++) {
    stepScene(alone, step);
  }

  const first = createScene();
  const second = createScene();
  for (let step = 0; step < 120; step++) {
    stepScene(first, step);
    stepScene(second, step + 1);
  }

  assert.deepEqual(readState(first.bodies), readState(alone.bodies));
});

test("StepFixed steps whole time steps and carries the remainder", () => {
  const world = createCorner();
  assert.equal(world.StepFixed(0.625, 0.25, 10, 10), 2);
  assert.equal(world.StepFixed(0.0625, 0.25, 10, 10), 0);
  assert.equal(world.StepFixed(0.0625, 0.25, 10, 10), 1);
  assert.equal(world.StepFixed(Number.NaN, 0.25, 10, 10), 0);
  assert.equal(world.StepFixed(-1, 0.25, 10, 10), 0);
  assert.equal(world.GetStepRemainder(), 0);
});

test("trigonometry matches fdlibm bit for bit", () => {
  assert.equal(b2Math.Sin(1), 0.8414709848078965);
  assert.equal(b2Math.Cos(1), 0.5403023058681398);
  assert.equal(b2Math.Sin(1e22), -0.8522008497671888);
  assert.equal(b2Math.Cos(1e22), 0.523214785395139);
  assert.equal(b2Math.Sin(-0), -0);
  assert.ok(Number.isNaN(b2Math.Sin(Infinity)));
  assert.equal(b2Math.Atan2(1, -1), 2.356194490192345);
  assert.equal(b2Math.Atan2(-0, -1), -Math.PI);
});
