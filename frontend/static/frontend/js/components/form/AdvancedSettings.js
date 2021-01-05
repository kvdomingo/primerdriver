import React, { Component } from "react";
import { MDBBtn as Button, MDBCollapse as Collapse } from "mdbreact";

export class AdvancedSettings extends Component {
  constructor(props) {
    super(props);
    this.state = {
      isOpen: false,
    };

    this.toggleCollapse = this.toggleCollapse.bind(this);
  }

  toggleCollapse() {
    this.setState(prevState => ({ isOpen: !prevState.isOpen }));
  }

  render() {
    return (
      <div className="accordion my-3" id="accordionAdvanced">
        <div className="card z-depth-0 bordered">
          <div className="card-header" id="headingAdvanced">
            <h5 className="mb-0 text-center">
              <Button
                className="btn btn-link"
                color="none"
                id="show-advanced"
                type="button"
                data-toggle="collapse"
                data-target="#advanced"
                aria-expanded="true"
                aria-controls="advanced"
                onClick={this.toggleCollapse}
              >
                Advanced settings
              </Button>
            </h5>
          </div>
          <Collapse
            id="advanced"
            aria-labelledby="headingAdvanced"
            data-parent="#accordionAdvanced"
            isOpen={this.state.isOpen}
          >
            <div className="card-body">
              <div className="row row-cols-3 text-center">
                <div className="col"></div>
                <div className="col">Min</div>
                <div className="col">Max</div>
              </div>
              <div className="row row-cols-3">
                <div className="col text-right">
                  Melting point (<sup>o</sup>C)
                </div>
                <div className="col">
                  <div className="md-outline">
                    <input
                      type="number"
                      id="Tm_range_min"
                      name="Tm_range_min"
                      className="form-control"
                      value={this.props.Tm_range_min}
                      onChange={this.props.handleChangeFloat}
                    />
                  </div>
                </div>
                <div className="col">
                  <div className="md-outline">
                    <input
                      type="number"
                      id="Tm_range_max"
                      name="Tm_range_max"
                      className="form-control"
                      value={this.props.Tm_range_max}
                      onChange={this.props.handleChangeFloat}
                    />
                  </div>
                </div>
              </div>
              <div className="row row-cols-3">
                <div className="col text-right">GC Content (%)</div>
                <div className="col">
                  <div className="md-outline">
                    <input
                      type="number"
                      id="gc_range_min"
                      name="gc_range_min"
                      className="form-control"
                      value={this.props.gc_range_min}
                      onChange={this.props.handleChangeFloat}
                    />
                  </div>
                </div>
                <div className="col">
                  <div className="md-outline">
                    <input
                      type="number"
                      id="gc_range_max"
                      name="gc_range_max"
                      className="form-control"
                      value={this.props.gc_range_max}
                      onChange={this.props.handleChangeFloat}
                    />
                  </div>
                </div>
              </div>
              <div className="row row-cols-3">
                <div className="col text-right">Length (bp)</div>
                <div className="col">
                  <div className="md-outline">
                    <input
                      type="number"
                      id="length_min"
                      name="length_min"
                      className="form-control"
                      value={this.props.length_min}
                      onChange={this.props.handleChangeInt}
                    />
                  </div>
                </div>
                <div className="col">
                  <div className="md-outline">
                    <input
                      type="number"
                      id="length_max"
                      name="length_max"
                      className="form-control"
                      value={this.props.length_max}
                      onChange={this.props.handleChangeInt}
                    />
                  </div>
                </div>
              </div>
              <div className="row row-cols-3">
                <div className="col text-right">5&apos; flanking region (bp)</div>
                <div className="col">
                  <div className="md-outline">
                    <input
                      type="number"
                      id="flank5_range_min"
                      name="flank5_range_min"
                      className="form-control"
                      value={this.props.flank5_range_min}
                      onChange={this.props.handleChangeInt}
                    />
                  </div>
                </div>
                <div className="col">
                  <div className="md-outline">
                    <input
                      type="number"
                      id="flank5_range_max"
                      name="flank5_range_max"
                      className="form-control"
                      value={this.props.flank5_range_max}
                      onChange={this.props.handleChangeInt}
                    />
                  </div>
                </div>
              </div>
              <div className="row">
                <div className="col text-right">3&apos; flanking region (bp)</div>
                <div className="col">
                  <div className="md-outline">
                    <input
                      type="number"
                      id="flank3_range_min"
                      name="flank3_range_min"
                      className="form-control"
                      value={this.props.flank3_range_min}
                      onChange={this.props.handleChangeInt}
                    />
                  </div>
                </div>
                <div className="col">
                  <div className="md-outline">
                    <input
                      type="number"
                      id="flank3_range_max"
                      name="flank3_range_max"
                      className="form-control"
                      value={this.props.flank3_range_max}
                      onChange={this.props.handleChangeInt}
                    />
                  </div>
                </div>
              </div>
              <div className="row row-cols-3">
                <div className="col text-right">5&apos; forward overlap (bp)</div>
                <div className="col">
                  <div className="md-outline">
                    <input
                      type="number"
                      id="forward_overlap5"
                      name="forward_overlap5"
                      className="form-control"
                      value={this.props.forward_overlap5}
                      onChange={this.props.handleChangeInt}
                    />
                  </div>
                </div>
                <div className="col"></div>
              </div>
              <div className="row row-cols-3">
                <div className="col text-right">3&apos; forward overlap (bp)</div>
                <div className="col">
                  <div className="md-outline">
                    <input
                      type="number"
                      id="forward_overlap3"
                      name="forward_overlap3"
                      className="form-control"
                      value={this.props.forward_overlap3}
                      onChange={this.props.handleChangeInt}
                    />
                  </div>
                </div>
                <div className="col"></div>
              </div>
              {this.props.mode === "PRO" ? (
                <div className="row row-cols-3">
                  <div className="col text-right">Expression system</div>
                  <div className="col">
                    <select
                      className="browser-default custom-select"
                      id="expression_system"
                      name="expression_system"
                      value={this.props.expression_system}
                      onChange={this.props.handleChange}
                      onMouseUp={e => {
                        e.target.name = "expression_system";
                        this.props.handleChange(e);
                      }}
                    >
                      {this.props.expressionList.map((exp, i) => (
                        <option value={exp} key={i}>
                          {exp}
                        </option>
                      ))}
                    </select>
                  </div>
                </div>
              ) : null}
              <div className="row row-cols-3">
                <div className="col text-right">Primer type</div>
                <div className="col">
                  <select
                    className="browser-default custom-select"
                    id="primer_mode"
                    name="primer_mode"
                    value={this.props.primer_mode}
                    onChange={this.props.handleChange}
                    onMouseUp={e => {
                      e.target.name = "primer_mode";
                      this.props.handleChange(e);
                    }}
                  >
                    <option value="complementary">Complementary</option>
                    <option value="overlapping">Overlapping</option>
                  </select>
                </div>
                <div className="col"></div>
              </div>
              <div className="custom-control custom-checkbox text-left mt-3">
                <input
                  type="checkbox"
                  className="custom-control-input"
                  id="terminate_gc"
                  name="terminate_gc"
                  checked={this.props.terminate_gc}
                  onChange={this.props.handleChangeChecked}
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
                  checked={this.props.center_mutation}
                  onChange={this.props.handleChangeChecked}
                />
                <label className="custom-control-label" htmlFor="center_mutation">
                  Mutation at center of primer
                </label>
              </div>
            </div>
          </Collapse>
        </div>
      </div>
    );
  }
}
