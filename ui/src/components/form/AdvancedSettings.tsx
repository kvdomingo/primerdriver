import {
  MDBCard as Card,
  MDBCardBody as CardBody,
  MDBCardHeader as CardHeader,
  MDBCol as Col,
  MDBCollapse as Collapse,
  MDBRow as Row,
} from "mdbreact";
import { useState } from "react";

export const defaultSettings = {
  Tm_range_min: 75,
  Tm_range_max: 85,
  gc_range_min: 40,
  gc_range_max: 60,
  length_min: 25,
  length_max: 45,
  flank5_range_min: 11,
  flank5_range_max: 21,
  flank3_range_min: 11,
  flank3_range_max: 21,
  forward_overlap5: 9,
  forward_overlap3: 9,
  terminate_gc: true,
  center_mutation: true,
  primer_mode: "complementary",
  expression_system: "Homo sapiens",
};

export type Settings = typeof defaultSettings;

const RANGES: [
  "Tm_range" | "gc_range" | "length" | "flank5_range" | "flank3_range",
  string,
][] = [
  ["Tm_range", "Melting point (°C)"],
  ["gc_range", "GC Content (%)"],
  ["length", "Length (bp)"],
  ["flank5_range", "5' flanking region (bp)"],
  ["flank3_range", "3' flanking region (bp)"],
];

interface AdvancedSettingsProps {
  value: Settings;
  onChange: (value: Settings) => void;
  mode?: string;
  expressionList?: string[];
}

function AdvancedSettings({
  value,
  onChange,
  mode,
  expressionList,
}: AdvancedSettingsProps) {
  const [isOpen, setIsOpen] = useState(false);

  function set(patch: Partial<Settings>) {
    onChange({ ...value, ...patch });
  }

  return (
    <div className="accordion my-3" id="accordionAdvanced">
      <Card className="bordered z-depth-0">
        <CardHeader id="headingAdvanced">
          <h5 className="mb-0 text-center">
            <button
              className="btn btn-link"
              color="none"
              id="show-advanced"
              type="button"
              aria-expanded={isOpen}
              aria-controls="advanced"
              onClick={() => setIsOpen(!isOpen)}
            >
              Advanced settings
            </button>
          </h5>
        </CardHeader>
        <Collapse id="advanced" aria-labelledby="headingAdvanced" isOpen={isOpen}>
          <CardBody>
            <Row className="row-cols-3 text-center">
              <Col />
              <Col>Min</Col>
              <Col>Max</Col>
            </Row>
            {RANGES.map(([prefix, label]) => {
              const minKey = `${prefix}_min` as const;
              const maxKey = `${prefix}_max` as const;
              const min = value[minKey];
              const max = value[maxKey];
              return (
                <Row className="row-cols-3" key={prefix}>
                  <Col className="text-right">{label}</Col>
                  <Col>
                    <div className="md-outline">
                      <input
                        type="number"
                        id={minKey}
                        name={minKey}
                        className="form-control"
                        value={min}
                        onChange={e => {
                          const v = Number.parseFloat(e.target.value);
                          set({ [minKey]: v < max ? v : max - 1 });
                        }}
                      />
                    </div>
                  </Col>
                  <Col>
                    <div className="md-outline">
                      <input
                        type="number"
                        id={maxKey}
                        name={maxKey}
                        className="form-control"
                        value={max}
                        onChange={e => {
                          const v = Number.parseFloat(e.target.value);
                          set({ [maxKey]: v > min ? v : min + 1 });
                        }}
                      />
                    </div>
                  </Col>
                </Row>
              );
            })}
            <hr />
            {(
              [
                ["forward_overlap5", "5' forward overlap (bp)"],
                ["forward_overlap3", "3' forward overlap (bp)"],
              ] as const
            ).map(([key, label]) => (
              <Row className="row-cols-3" key={key}>
                <Col className="text-right">{label}</Col>
                <Col>
                  <div className="md-outline">
                    <input
                      type="number"
                      id={key}
                      name={key}
                      className="form-control"
                      value={value[key]}
                      onChange={e => set({ [key]: Number.parseInt(e.target.value) })}
                    />
                  </div>
                </Col>
                <Col />
              </Row>
            ))}
            {mode === "PRO" ? (
              <Row className="row-cols-3">
                <Col className="text-right">Expression system</Col>
                <Col>
                  <select
                    className="browser-default custom-select"
                    id="expression_system"
                    name="expression_system"
                    value={value.expression_system}
                    onChange={e => set({ expression_system: e.target.value })}
                  >
                    {(expressionList ?? []).map(exp => (
                      <option value={exp} key={exp}>
                        {exp}
                      </option>
                    ))}
                  </select>
                </Col>
              </Row>
            ) : null}
            <Row className="row-cols-3">
              <Col className="text-right">Primer type</Col>
              <Col>
                <select
                  className="browser-default custom-select"
                  id="primer_mode"
                  name="primer_mode"
                  value={value.primer_mode}
                  onChange={e => set({ primer_mode: e.target.value })}
                >
                  <option value="complementary">Complementary</option>
                  <option value="overlapping">Overlapping</option>
                </select>
              </Col>
              <Col />
            </Row>
            <Row className="row-cols-3">
              <Col />
              <Col>
                <div className="custom-control custom-checkbox mt-3 text-left">
                  <input
                    type="checkbox"
                    className="custom-control-input"
                    id="terminate_gc"
                    name="terminate_gc"
                    checked={value.terminate_gc}
                    onChange={e => set({ terminate_gc: e.target.checked })}
                  />
                  <label className="custom-control-label" htmlFor="terminate_gc">
                    Terminates in G/C
                  </label>
                </div>
                <div className="custom-control custom-checkbox text-left">
                  <input
                    type="checkbox"
                    className="custom-control-input"
                    id="center_mutation"
                    name="center_mutation"
                    checked={value.center_mutation}
                    onChange={e => set({ center_mutation: e.target.checked })}
                  />
                  <label className="custom-control-label" htmlFor="center_mutation">
                    Mutation at center of primer
                  </label>
                </div>
              </Col>
              <Col />
            </Row>
          </CardBody>
        </Collapse>
      </Card>
    </div>
  );
}

export { AdvancedSettings };
