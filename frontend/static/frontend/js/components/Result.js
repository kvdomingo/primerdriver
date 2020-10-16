import React, { Component } from "react";
import { Link } from "react-router-dom";
import {
  MDBTable as Table,
  MDBTableHead as TableHead,
  MDBTableBody as TableBody,
  MDBIcon as Icon,
  MDBTypography as Typography,
} from "mdbreact";

export default class Result extends Component {
  header() {
    return (
      <div>
        <Link to="/" className="btn btn-blue-grey mb-4 mr-3" id="back">
          <Icon fas icon="arrow-left" className="mr-3" />
          main menu
        </Link>
        <Typography tag="h2" className="d-md-inline ml-md-3">
          {(() => {
            if (this.props.mode === "CHAR") return "Primer characterization";
            else if (this.props.mode === "DNA") return "DNA-based primer design";
            else if (this.props.mode === "PRO") return "Protein-based primer design";
            else return null;
          })()}
        </Typography>
      </div>
    );
  }

  render() {
    if (typeof this.props.results === "object") {
      if (this.props.mode === "CHAR")
        return (
          <div>
            {this.header()}
            <Table responsive bordered size="sm" hover className="text-nowrap">
              <TableBody>
                {Object.keys(this.props.results).map((key, i) => (
                  <tr key={i}>
                    <th scope="row">{key}</th>
                    <td>{this.props.results[key]["1"]}</td>
                  </tr>
                ))}
              </TableBody>
            </Table>
          </div>
        );
      else
        return (
          <div>
            {this.header()}
            <Typography tag="h2" variant="h2-responsive" className="mx-md-2 my-4">
              {Object.keys(this.props.results).length} results
            </Typography>
            {Object.keys(this.props.results).map((item, i) => (
              <Table responsive className="text-nowrap" size="sm" bordered hover key={i}>
                <TableHead>
                  <tr>
                    <th scope="col" />
                    <th scope="col">{`Primer ${i + 1}`}</th>
                  </tr>
                </TableHead>
                <TableBody>
                  {Object.keys(this.props.results[(i + 1).toString()]).map((key, j) => (
                    <tr key={j}>
                      <th scope="row">{key}</th>
                      <td>{this.props.results[(i + 1).toString()][key]}</td>
                    </tr>
                  ))}
                </TableBody>
              </Table>
            ))}
          </div>
        );
    } else {
      return (
        <div className="text-center">
          {this.header()}
          <p>
            Oops! Something went wrong on the server. Please try again later, or{" "}
            <a href="https://github.com/kvdomingo/primerdriver/issues" target="_blank" rel="noopener noreferrer">
              report an issue
            </a>
            .
          </p>
        </div>
      );
    }
  }
}
